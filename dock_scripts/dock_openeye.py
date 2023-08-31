#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 17:17:05 2021

@author: alexehaywood
"""

# Import General Proccess Modules
import argparse
import os
import sys
import time
import re
import operator
import subprocess as sp
from collections import defaultdict

# Import Openeye Modules
from openeye import oechem
from openeye import oeomega
from openeye import oequacpac
from openeye import oedocking

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import AllChem

NUM_POSES = 10000

# Find Openeye licence
try:
    print('license file found: %s' % os.environ['OE_LICENSE'])
except KeyError:
    print('license file not found, please set $OE_LICENSE')
    sys.exit('Critical Error: license not found')

# Find obabel and tautomers locations
obabel_loc = sp.run(['which', 'obabel'], stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True).stdout.strip()
oetautomers_loc = sp.run(['which', 'tautomers'], stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True).stdout.strip()

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
    # AllChem.EmbedMolecule(enhanced_mol_block, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
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
        else:
            # OESetNeutralpHModel should return True even if ligand is unchanged
            sys.exit('ERROR: OESetNeutralpHModel returned False for ligand: {}'.format(ligand_prefix))

    ifs.close()
    ofs.close()
    #os.rename(ligand_sdf_fpath[:-4]+'_pH74.sdf', ligand_sdf_fpath)
    
    
    sp.call('obabel ' + ligand_sdf_fpath + ' -O ' + ligand_prefix + '.svg', shell=True)
    
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
    
    # conformer options
    omega_opts = oeomega.OEOmegaOptions()
    omega_opts.SetEnergyWindow(50)
    omega_opts.SetMaxSearchTime(600)
    omega_opts.SetSearchForceField(7)
    omega_opts.SetRMSThreshold(0.25) 
    omega_opts.SetMaxConfs(10000) 
    omega_opts.SetStrictStereo(True) # change to True if stereocenters are known
    
    omega = oeomega.OEOmega(omega_opts)
    
    # input molecule file
    ifs = oechem.oemolistream()
    ifs.open(ligand_sdf_fpath)
    
    # output molecule file
    ofs = oechem.oemolostream()
    ofs.open(ligand_oeb_fpath); print('file %s generated' %ligand_oeb_fpath)
    
    # tautomer options
    #tautomer_opts = oequacpac.OETautomerOptions()
    #tautomer_opts.SetMaxSearchTime(300)
    #tautomer_opts.SetRankTautomers(True)
    #tautomer_opts.SetCarbonHybridization(False)
    # pKa_norm = True
    
    # enantiomer options
    #flipper_opts = oeomega.OEFlipperOptions()
    #flipper_opts.SetMaxCenters(12)
    # Preserve any defined stereocentres:
    #flipper_opts.SetEnumSpecifiedStereo(False)
    #flipper_opts.SetEnumNitrogen(True)
    #flipper_opts.SetWarts(False)
    
    # generate tautomers, enantiomers and conformers
    for mol in ifs.GetOEMols():
        
        #for tautomer in oequacpac.OEEnumerateTautomers(mol, tautomer_opts):
        # for tautomer in oequacpac.OEGetReasonableTautomers(mol, tautomer_opts, pKa_norm):
            
            # comment out the next line if specific stereocenter is know and encoded in the original mol file
            #for enantiomer in oeomega.OEFlipper(tautomer, flipper_opts):
                #ligand = oechem.OEMol(enantiomer)
                # Add ligand name for sdf file
        mol.SetTitle(ligand_prefix)
        
        # ligand_check = oechem.OEGraphMol()
        # oechem.OESmilesToMol(ligand_check, oechem.OEMolToSmiles(ligand))
        # parent_search = oechem.OESubSearch('C') 
        # oechem.OEPrepareSearch(ligand_check, parent_search)
        
        #if nonallowed_fragment_tautomers(mol):
        ret_code = omega.Build(mol)
            
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            oechem.OEWriteMolecule(ofs, mol)

                        
    ifs.close()
    ofs.close()
    

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
    #ifs.open(ligand_oeb_fpath); print('file %s opened' %ligand_oeb_fpath); os.system('ls %s' %ligand_oeb_fpath)
    ifs.open(ligand_oeb_fpath)

    # receptor file
    rfs = oechem.oeifstream()
    rfs.open(receptor_oedu_fpath)

    # output molecule file
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
    
    # output file to save docking scoress to
    output_file = open(output_fpath, 'w+')
    
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
    output_file.write('Docking Score: ')
    for ligand in ifs.GetOEMols():
        print("docking", ligand.GetTitle())
        # docked_mol = oechem.OEGraphMol()
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
    
    if oechem.OEMolBase.GetEnergy(best_pose) != best_docking_score:
        print('best pose energy = ', oechem.OEMolBase.GetEnergy(best_pose) )
        print('best docking score = ', best_docking_score)
        sys.exit('ERROR: Check best pose.')
    oechem.OEWriteConstMolecule(ofs_sdf_best, best_pose)
           
    ifs.close()
    ofs.close()
    ofs_sdf.close()
    ofs_sdf_best.close()
    rfs.close()
    output_file.close()
        
    return best_docking_score


def oe_dock(ligand_mol_fpath, receptor_oedu_fpath):
    """
    Prepare the ligand, generate tautomers/enantiomers/conformers and dock into
    the receptor(s).

    Parameters
    ----------
    ligand_mol_fpath : str
        ligand file path (.mol or .sdf).
    receptor_oeb_fpaths : str
        receptor(s) .oedu file path(s).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    # This is the function that gets invoked from C++
    start_t = time.time()
    
    #receptors = re.split('(?<!auto)[+\-,](?!docker)', receptor_oedu_fpaths)
    
    ligand_prefix = ligand_mol_fpath[:-4]
    ligand_oeb_fpath = ligand_prefix + '.oeb'
    
    #molecule_prep(ligand_mol_fpath, ligand_prefix)
    
    oe_conformer_generation(ligand_prefix)	
    
    dscore_dict = defaultdict()
    #for receptor_oedu_fpath in receptors:
    receptor_prefix = receptor_oedu_fpath[:-5]
    docked_oeb_fpath = ligand_prefix + '_' + receptor_prefix + '-docked.oeb'
    output_fpath = '{}_{}.score'.format(ligand_prefix, receptor_prefix)
    
    best_docking_score = dock(ligand_oeb_fpath, receptor_oedu_fpath, docked_oeb_fpath, output_fpath)
    
    dscore_dict[receptor_prefix] = best_docking_score
    
    # Write the properties in the *receptor.oe.score file
    output_file = open(output_fpath, 'a+')
    output_file.write(
        '\nTime: {} \nBest_Docking_Score: {}'.format(
            time.time() - start_t, best_docking_score)
        )
    output_file.close()
    
    #os.remove(ligand_oeb_fpath)

    return best_docking_score


def convert_pdb_to_oeb(input, output):

    ifs = oechem.oemolistream(input)
    ofs = oechem.oemolostream(output)

    for mol in ifs.GetOEGraphMols():
        oechem.OEWriteMolecule(ofs, mol)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Dock ligand to receptor using Openeye')

    parser.add_argument('-l','--ligand', help='Ligand file', required=True)
    parser.add_argument('-r','--receptor', help='Receptor file', required=True)

    args = parser.parse_args()

    ligand_mol_fpath = args.ligand
    receptor_oedu_fpath = args.receptor

    oe_dock(ligand_mol_fpath, receptor_oedu_fpath)
