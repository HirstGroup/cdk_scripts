#!/usr/bin/python

###!/cm/shared/apps/python/2.7.9/bin/python
####!/usr/bin/python

############################################################################################################################################
# The MIT License (MIT)                                                                                                                    #
#                                                                                                                                          #
# Copyright (c) 2017 Steven Oatley & Dino Oglic                                                                                            #
#                                                                                                                                          #
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files         #
# (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge,      #
# publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,   #
# subject to the following conditions:                                                                                                     #
#                                                                                                                                          #
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.           #
#                                                                                                                                          #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF       #
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR  #
# ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH   #
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                                                               #
############################################################################################################################################
# IMPORTANT NOTE: The code below relies on the code from a python library that is subject to a different license. Please consult the       #
# corresponding license before using the Software.                                                                                         #
############################################################################################################################################
# This code was mostly written by Steve Oatley, with a few modifications by Arnaldo F. Silva Filho (AFS), arsfilho@gmail.com               #
############################################################################################################################################


# =============================================================================
# Import Modules
# =============================================================================

# Import General Proccess Modules
import os
import sys
import time
import re
import operator
import subprocess as sp
import math
import pandas as pd
from collections import defaultdict

# Import Openeye Modules
from openeye import oechem, oeomega, oedocking, oequacpac

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import AllChem

# Add path so the predictive_models and properties modules can be found
head_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(head_path)

# Add to sys.argv to solve tensorflow loading error due to calling code from C++
if not hasattr(sys, "argv") or not sys.argv:
    sys.argv = [""]

# Import QSAR Models
from predictive_models.ml_model_gcnn_ens import Ensemble_Model_DC

# Import properties modules
from properties import oeMolProp, num_atomatic_rings, num_chiral_centres, \
    num_lipinsky_donors, num_lipinsky_acceptors, molecular_weight, num_rot_bond


# =============================================================================
# Define Thresholds
# =============================================================================

#   NUM_POSES:          Maximum number of top scoring docked poses to return
NUM_POSES = 15

# 	CHIRAL_THRESHOLD:   Maximum number of Chiral Centers (<= 2)
CHIRAL_THRESHOLD = 2
#	PSA_TRESHOLD:       Polar surface area (<= 140)
PSA_TRESHOLD = 140
#	PFI_TRESHOLD:       Property Forecast Index (LOGP + # of aromatic rings < 8)
PFI_TRESHOLD = 8
#   ROTBOND_THRESHOLD:  Number of rotatable bonds (<= 7)
ROTBOND_THRESHOLD = 7	
#   WEIGHT_TRESHOLD:    Maximum MW in Daltons (<= 500)
WEIGHT_TRESHOLD = 500
#   H_DON_TRESHOLD:     Maximum number of hydrogen donors (<= 5)
H_DON_TRESHOLD = 5
#   H_ACC_TRESHOLD:     Maximum number of hydrogen acceptors (<= 10)
H_ACC_TRESHOLD = 10
#   LOGP_TRESHOLD:      Water/Octanol Partition Coefficient (0.5-5.0)
LOGP_TRESHOLD_UP = 5
LOGP_TRESHOLD_LOW = 0.5
#   DOCKING_THRESHOLD:  Fred docking score (<= -6)
DOCKING_THRESHOLD = -6


# =============================================================================
# Initialisation
# =============================================================================

# Find Openeye licence
try:
    print('license file found: %s' % os.environ['OE_LICENSE'])
except KeyError:
    print('license file not found, please set $OE_LICENSE')
    sys.exit('Critical Error: license not found')

# Find obabel and tautomers locations
obabel_loc = sp.run(
    ['which', 'obabel'], stdout=sp.PIPE, stderr=sp.PIPE,
    universal_newlines=True).stdout.strip()
oetautomers_loc = sp.run(
    ['which', 'tautomers'], stdout=sp.PIPE, stderr=sp.PIPE,
    universal_newlines=True).stdout.strip()

# AFS: Define Pan-assay interference compounds (PAINS), read from the PAINS.csv
#      file in the directory
# https://pubs.acs.org/doi/10.1021/jm5019093
pains_df = pd.read_csv(os.path.join(head_path, 'PAINS.csv'))
pains_df['mol'] = pains_df['Smarts'].apply(Chem.MolFromSmarts)
pains_df['smiles'] = pains_df['mol'].apply(Chem.MolToSmiles)
PAINS_fragment_list = pains_df['smiles']

# Load the QSAR models
model_pIC50 = Ensemble_Model_DC(os.path.join(head_path, 'pIC50.pk'))


# =============================================================================
# Define Modules
# =============================================================================

# AFS: Filter out PAINS compunds, the PAINS list is in this DIR, as PAINS.csv
def PAINS_filter(molecule):
    for fragment in PAINS_fragment_list:
        fragment_search = oechem.OESubSearch(fragment)
        oechem.OEPrepareSearch(molecule, fragment_search)
        if fragment_search.SingleMatch(molecule):
            return False
            break

    return True

# AFS: Filter out undesirable tautomers
def nonallowed_fragment_tautomers(molecule):

    fragment_list = ["N=CO",
                     "C=NCO",
                     "C(O)=NC",
                     "c=N",
                     "N=c",
                     "C=n",
                     "n=C"]

    for fragment in fragment_list:
        fragment_search = oechem.OESubSearch(fragment)
        oechem.OEPrepareSearch(molecule, fragment_search)
        if fragment_search.SingleMatch(molecule):
            return False
            break

    return True

def rdkit_mol_enhance(fpath, output_fpath):

    mol_block = Chem.MolFromMolFile(fpath)
    enhanced_mol_block = Chem.AddHs(mol_block)
    # IF rdkit version > 2013.09.1 THEN
    # AllChem.EmbedMolecule(
    #   enhanced_mol_block, useExpTorsionAnglePrefs=True,
    #   useBasicKnowledge=True)
    # ELSE
    AllChem.EmbedMolecule(enhanced_mol_block)
    enhanced_mol_block = Chem.MolToMolBlock(enhanced_mol_block)

    sdf = open(output_fpath, 'w')
    sdf.write(enhanced_mol_block)
    sdf.close()

def molecule_prep(ligand_mol_fpath, ligand_prefix):
    """
    Prepare molecule by generating 3D coordinates, calculating ionisation state
    at pH 7.4 and saving an image of the molecule.
    Parameters
    ----------
    ligand_mol_fpath : str
        ligand file path (.mol or .sdf).
    ligand_prefix : str
        ligand file name without file extension.
    """
    ligand_sdf_fpath = '{}.sdf'.format(ligand_prefix)

    rdkit_mol_enhance(ligand_mol_fpath, ligand_sdf_fpath)

    ifs = oechem.oemolistream()
    ifs.open(ligand_sdf_fpath)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_SDF)
    ofs.open(ligand_sdf_fpath[:-4]+'_pH74.sdf')
    
    mol = oechem.OEGraphMol()
    while oechem.OEReadMolecule(ifs, mol):
        if oequacpac.OESetNeutralpHModel(mol):
            # Make any additional hydrogens added by OESetNeutralpHModel 
            # explicit and add 2D coordinates for the sdf file
            oechem.OEAddExplicitHydrogens(mol)
            oechem.OESet2DHydrogenGeom(mol)
            oechem.OEWriteMolecule(ofs, mol)

    ifs.close()
    ofs.close()
    os.rename(ligand_sdf_fpath[:-4]+'_pH74.sdf', ligand_sdf_fpath)

    sp.call('obabel ' + ligand_sdf_fpath + ' -O ' + ligand_prefix + '.svg',
            shell=True)

    return 0

def oe_conformer_generation(ligand_prefix):
    """
    Generate tautomers, enantiomers and conformers of the ligand. Ensure ligand
    sdf file is in directory.
    Parameters
    ----------
    ligand_prefix : str
        ligand file name without file extension.
    Returns
    -------
    filtered_by_pains : bool
        If True, the ligand was filtered out because it contains a PAINS fragment.
        If False, the ligand was not filtered out by the PAINS filter.
    """
    ligand_sdf_fpath = '{}.sdf'.format(ligand_prefix)
    ligand_oeb_fpath = '{}.oeb'.format(ligand_prefix)

    # input molecule file
    ifs = oechem.oemolistream()
    ifs.open(ligand_sdf_fpath)

    # output molecule file
    ofs = oechem.oemolostream()
    ofs.open(ligand_oeb_fpath)

    # conformer options
    omega_opts = oeomega.OEOmegaOptions()
    omega_opts.SetEnergyWindow(20)
    omega_opts.SetMaxSearchTime(600)
    omega_opts.SetSearchForceField(7)
    omega_opts.SetRMSThreshold(0.5)
    omega_opts.SetMaxConfs(1000)
    omega_opts.SetStrictStereo(False)

    omega = oeomega.OEOmega(omega_opts)

    # tautomer options
    tautomer_opts = oequacpac.OETautomerOptions()
    tautomer_opts.SetMaxSearchTime(300)
    tautomer_opts.SetRankTautomers(True)
    tautomer_opts.SetCarbonHybridization(False)
    # pKa_norm = True

    # enantiomer options
    flipper_opts = oeomega.OEFlipperOptions()
    flipper_opts.SetMaxCenters(12)
    flipper_opts.SetEnumSpecifiedStereo(True)
    flipper_opts.SetEnumNitrogen(True)
    flipper_opts.SetWarts(False)

    filtered_by_pains = []
    # generate tautomers, enantiomers and conformers
    for mol in ifs.GetOEMols():

        for tautomer in oequacpac.OEEnumerateTautomers(mol, tautomer_opts):
        # oequacpac.OEGetReasonableTautomers(mol, tautomer_opts, pKa_norm)
        # doesn't work properly

            for enantiomer in oeomega.OEFlipper(tautomer, flipper_opts):
                ligand = oechem.OEMol(enantiomer)

                if nonallowed_fragment_tautomers(ligand) and PAINS_filter(ligand):
                    ret_code = omega.Build(ligand)

                    if ret_code == oeomega.OEOmegaReturnCode_Success:
                        oechem.OEWriteMolecule(ofs, ligand)

                if PAINS_filter(ligand):
                    filtered_by_pains.append(False)
                else:
                    filtered_by_pains.append(True)

    ifs.close()
    ofs.close()

    if any(i is False for i in filtered_by_pains):
        return False
    else:
        return True

def dock(ligand_oeb_fpath, receptor_oedu_fpath, docked_oeb_fpath, output_fpath):
    """
    Dock ligand into receptor.
    Parameters
    ----------
    ligand_oeb_fpath : str
        ligand .oeb file path.
    receptor_oedu_fpath : str
        receptor .oedu file path.
    docked_oeb_fpath : str
        docked molecule .oeb file path.
    output_fpath : str
        file path to save docking results.
    Returns
    -------
    best_docking_score : float
        best docking score out of all tautomers, enantiomers and conformers in
        the ligand_oeb_fpath file.
    """
    # input molecule file
    ifs = oechem.oemolistream()
    ifs.open(ligand_oeb_fpath)

    # receptor file
    rfs = oechem.oeifstream()
    rfs.open(receptor_oedu_fpath)

    # output molecule oeb file
    ofs = oechem.oemolostream()
    ofs.open(docked_oeb_fpath)

    # output molecule sdf file
    ofs_sdf = oechem.oemolostream()
    ofs_sdf.SetFormat(oechem.OEFormat_SDF)
    ofs_sdf.open(docked_oeb_fpath[:-4]+'.sdf')

    # output molecule best pose sdf file for convenience
    ofs_sdf_best = oechem.oemolostream()
    ofs_sdf_best.SetFormat(oechem.OEFormat_SDF)
    ofs_sdf_best.open(docked_oeb_fpath[:-4]+'_best_pose.sdf')

    # output file to save docking scores to
    output_file = open(output_fpath, 'a')

    # read receptor
    receptor = oechem.OEDesignUnit()
    oechem.OEReadDesignUnit(rfs, receptor)
    receptor.HasReceptor()

    # docking parameters
    dock_method = oedocking.OEDockMethod_Chemgauss4
    dock_resolution = oedocking.OESearchResolution_High

    # initialise docking
    dock = oedocking.OEDock(dock_method, dock_resolution)
    dock.Initialize(receptor)

    # dock molecules
    best_docking_score = 1000000000
    for ligand in ifs.GetOEMols():
        print("docking", ligand.GetTitle())

        docked_mol = oechem.OEMol()
        dock.DockMultiConformerMolecule(docked_mol, ligand, NUM_POSES)

        sdtag = oedocking.OEDockMethodGetName(dock_method)
        oedocking.OESetSDScore(docked_mol, dock, sdtag)
        dock.AnnotatePose(docked_mol)
        oechem.OEWriteMolecule(ofs, docked_mol)

        docking_score = oechem.OEMolBase.GetEnergy(docked_mol)
        output_file.write('{} '.format(docking_score))

        if docking_score < best_docking_score:
            best_docking_score = docking_score
            best_pose = oechem.OEMol()
            oechem.OECopyMol(best_pose, docked_mol.GetConfs().next())

        oechem.OEWriteConstMolecule(ofs_sdf, docked_mol)

    output_file.write('\n')

    oechem.OEWriteConstMolecule(ofs_sdf_best, best_pose)

    ifs.close()
    ofs.close()
    ofs_sdf.close()
    ofs_sdf_best.close()
    rfs.close()
    output_file.close()
        
    return best_docking_score

def oe_dock(ligand_mol_fpath, receptor_oedu_fpaths):
    """
    Prepare the ligand, generate tautomers/enantiomers/conformers and dock into
    the receptor(s).
    Parameters
    ----------
    ligand_mol_fpath : str
        ligand file path (.mol or .sdf).
    receptor_oedu_fpaths : str
        receptor(s) .oedu file path(s).
    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    # This is the function that gets invoked from C++
    start_t = time.time()

    receptors = re.split('(?<!auto)[+\-,](?!docker)', receptor_oedu_fpaths)

    ligand_prefix = os.path.splitext(ligand_mol_fpath)[0]
    ligand_sdf_fpath = ligand_prefix + '.sdf'
    ligand_oeb_fpath = ligand_prefix + '.oeb'

    output_fpath = '{}.score'.format(ligand_prefix)

    # Generate 3D coordinates, save to sdf file, calculate ionisation state at
    # pH 7.4 and save image
    molecule_prep(ligand_mol_fpath, ligand_prefix)

    # Perform QSAR predictions and calculate synthetic accessibility scores
    ifs = oechem.oemolistream()
    ifs.open(ligand_sdf_fpath)
    for mol in ifs.GetOEMols():
        ligand_smi = oechem.OEMolToSmiles(mol)

        predicted_pIC50, warnings_pIC50 = model_pIC50.predict(ligand_smi)

    # Pre-compute logP, PSA, PFI, H_donors, H_acceptors, molecular weight,
    # number of rot bonds to pre-filter compounds
    logp, PSA = oeMolProp(ligand_sdf_fpath)
    n_aromatic_rings = num_atomatic_rings(ligand_smi)
    PFI = n_aromatic_rings + logp
    n_chiral = num_chiral_centres(ligand_sdf_fpath)
    h_don = num_lipinsky_donors(ligand_sdf_fpath)
    h_acc = num_lipinsky_acceptors(ligand_sdf_fpath)
    mol_weight = molecular_weight(ligand_sdf_fpath)
    n_rot_bonds = num_rot_bond(ligand_sdf_fpath)

    # Write properties to output file
    output_file = open(output_fpath, 'w+')
    output_file.write(
        'Ligand: {}\n'.format(os.path.basename(ligand_prefix))
        + 'Predicted_pIC50: {}\n'.format(predicted_pIC50)
        + 'logP: {}\n'.format(logp)
        + 'N_Aromatic_Rings: {}\n'.format(n_aromatic_rings)
        + 'PFI: {}\n'.format(PFI)
        + 'PSA: {}\n'.format(PSA)
        + 'Molecular_Weight: {}\n'.format(mol_weight)
        + 'N_Chiral_Centres: {}\n'.format(n_chiral)
        + 'N_Acceptors: {}\n'.format(h_acc)
        + 'N_Donors: {}\n'.format(h_don)
        + 'N_Rotatable_Bonds: {}\n'.format(n_rot_bonds)
        )
    output_file.close()

    # If ligand is filterd out, just return a very high score
    if logp > LOGP_TRESHOLD_UP or logp < LOGP_TRESHOLD_LOW \
            or n_chiral > 2 \
            or PSA > PSA_TRESHOLD or PFI > PFI_TRESHOLD \
            or mol_weight > WEIGHT_TRESHOLD \
            or h_don > H_DON_TRESHOLD or h_acc > H_ACC_TRESHOLD \
            or n_rot_bonds > ROTBOND_THRESHOLD:

        # Write filtered out by filters to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            'Filtered_Out_by_Filters: True\n'
            + 'Time: {}\n'.format(time.time() - start_t)
            )
        output_file.close()

        # Return high MPO score
        return 1000000000.0

    else:
        # Write filtered out by filters to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            'Filtered_Out_by_Filters: False\n'
            )
        output_file.close()

    # Generate conformers, tautomers and enantiomers
    filtered_by_pains = oe_conformer_generation(ligand_prefix)
    
    # If ligand if filtered out by PAINS, return a very high score
    if filtered_by_pains is True:
        # Write filtered out by PAINS to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            'Filtered_Out_by_Pains: True\n'
            + 'Time: {}\n'.format(time.time() - start_t)
            )
        output_file.close()

        # Return high MPO score
        return 1000000000.0

    else:
        # Write filtered out by PAINS to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            'Filtered_Out_by_Pains: False\n'
            )
        output_file.close()

    # Calculate docking score and MPO function
    score_dict = defaultdict()
    for receptor_oedu_fpath in receptors:
        receptor_prefix = os.path.basename(os.path.splitext(receptor_oedu_fpath)[0])
        docked_oeb_fpath = ligand_prefix + '_' + receptor_prefix + '-docked.oeb'
        # output_fpath = '{}_{}.score'.format(ligand_prefix, receptor_prefix)
    
        # Write the receptor name to output file
        output_file = open(output_fpath, 'a')
        output_file.write('\nReceptor: {}\n'.format(receptor_prefix))
        output_file.write('{}_Docking_Scores: '.format(receptor_prefix))
        output_file.close()

        best_docking_score = dock(
            ligand_oeb_fpath, receptor_oedu_fpath, docked_oeb_fpath,
            output_fpath
            )

        if best_docking_score > DOCKING_THRESHOLD:
            # Write filtered out by docking to output file
            output_file = open(output_fpath, 'a')
            output_file.write(
                '{}_Best_Docking_Score: {}\n'.format(receptor_prefix, best_docking_score)
                + '{}_Filtered_Out_by_Docking: True\n'.format(receptor_prefix)
                + 'Time: {}\n'.format(time.time() - start_t)
                )
            output_file.close()

            # Return high MPO score
            return 1000000000.0
        
        else:
            # MPO equation: First term potency, second term sigmoidal equation
            # for solubility. Add more parameters below.
            mpo = (- predicted_pIC50) * ( 1 / (1 + math.exp(PFI - 8)) )
    
            score_dict[receptor_prefix] = mpo
    
            # Write the best docking score to output file
            output_file = open(output_fpath, 'a')
            output_file.write(
                '{}_Best_Docking_Score: {}\n'.format(receptor_prefix, best_docking_score)
                + '{}_Filtered_Out_by_Docking: False\n\n'.format(receptor_prefix)
                )
            output_file.close()

    os.remove(ligand_oeb_fpath)

    # If more than one receptor
    if len(receptors) > 1:
        ops = {'+': operator.add, ',': operator.add, '-': operator.sub}
        delta = 0
        recep_separators = re.findall(
            '(?<!auto)[+\-,](?!docker)', receptor_oedu_fpaths
            )
        recep_separators.insert(0, '+')
        receptor_prefix = receptor_oedu_fpath.split('/')[-1][:-4]
        print([receptor.split('/')[-1][:-4] for receptor in receptors])
        scores = [score_dict[receptor.split('/')[-1][:-4]] for receptor in receptors]
        delta = [ops[recep_sep](delta, score) for score, recep_sep in zip(scores, recep_separators)]

        # Write the MPO score to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            + 'MPO: {}\n'.format(delta)
            + 'Time: {}\n'.format(time.time() - start_t)
            )
        output_file.close()
    
        return mpo

    # If one receptor return MPO score
    else:
        # Write the MPO score to output file
        output_file = open(output_fpath, 'a')
        output_file.write(
            'MPO: {}\n'.format(mpo)
            + 'Time: {}\n'.format(time.time() - start_t)
            )
        output_file.close()
        
        return mpo


def setup_test_vars():
    
    ligand_mol_fpath = '/home/alexehaywood/Documents/active_search/oracle/testing/phenylisoxazole.mol'
    receptor_oedu_fpaths = '/home/alexehaywood/Documents/active_search/docking/docking_code/receptors/receptor_contraints.oedu'
    
    return ligand_mol_fpath, receptor_oedu_fpaths


if __name__ == '__main__':
    
    # ligand_mol_fpath, receptor_oedu_fpaths = setup_test_vars()
    # oe_dock(ligand_mol_fpath, receptor_oedu_fpaths)

    oe_dock(sys.argv[1], sys.argv[2])