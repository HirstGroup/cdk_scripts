#!/usr/bin/env python3

# Code to generate molecule sdf files for GNINA docking.

# Based on dock_molecule.py from Alexe, but with docking 
# functions removed and option to turn TorsionDrive off, 
# ready for docking with GNINA, which will search torsion 
# angles itself.  Should also preserve any defined 
# stereochemistry while enumerating over undefined 
# stereochemistry.

# Import General Proccess Modules
import os
import sys
import time
import re
import operator
import subprocess as sp
from collections import defaultdict
import argparse

# Import Openeye Modules
from openeye import oechem
from openeye import oeomega
from openeye import oequacpac

# Import RDKit Tools
from rdkit import Chem
from rdkit.Chem import AllChem

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

def rdkit_mol_enhance(ligand_smiles, output_fpath):
    
    # mol_block = Chem.MolFromMolFile(fpath)
    mol_block = Chem.MolFromSmiles(ligand_smiles)
    enhanced_mol_block = Chem.AddHs(mol_block)
    # IF rdkit version > 2013.09.1 THEN
    # AllChem.EmbedMolecule(enhanced_mol_block, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    # ELSE

    # Don't embed molecule here, otherwise it causes undefined stereochemistry to become defined:
    # AllChem.EmbedMolecule(enhanced_mol_block)

    enhanced_mol_block = Chem.MolToMolBlock(enhanced_mol_block)
    sdf = open(output_fpath, 'w')
    sdf.write(enhanced_mol_block)
    sdf.close()

def molecule_prep(ligand_smiles, ligand_prefix, ph=True):
    """
    Prepare molecule by generating 3D coordinates, calculating ionisation state 
    at pH 7.4 and saving an image of the molecule.

    Parameters
    ----------
    ligand_mol_fpath : str
        ligand file path (.mol or .sdf)
    ligand_prefix : str
        ligand file name without file extension
    ph : bool
        convert ligand to protonation state at pH 7.4
    """
    ligand_sdf_fpath = '{}.sdf'.format(ligand_prefix)

    rdkit_mol_enhance(ligand_smiles, ligand_sdf_fpath)   

    if ph:
        print('Convert to pH 7.4 (OpenEye): {}'.format(ph))
        ifs = oechem.oemolistream()                           
        ifs.open(ligand_sdf_fpath) 
    
        ofs = oechem.oemolostream()                            
        ofs.SetFormat(oechem.OEFormat_SDF)    

        # Need to save ligand in state at pH 7.4 to file:
        ofs.open(ligand_sdf_fpath[:-4]+'_pH74.sdf')
    
        mol = oechem.OEGraphMol()
        while oechem.OEReadMolecule(ifs, mol):
            # OESetNeutralpHModel returns True if it runs successfully, even if mol is unchanged
            if oequacpac.OESetNeutralpHModel(mol):
                # Need to make any hydrogens added during pH adjustment 
                # explicit so that they will be added to the sdf file
                oechem.OEAddExplicitHydrogens(mol) #, set3D=False)
                
                # Explicit hydrogens will be given same coordinates as 
                # associated heavy atom, need to generate new coordinates
                # This seems to reset stereochemistry so should not be
                # used here if molecules have some stereochemistry.  Not
                # clear why this is, but possibly it means that OpenEye
                # considers the molecule as having proper coordinates 
                # which then overrides any stereochemistry flags in the
                # sdf file:
                #oechem.OESet2DHydrogenGeom(mol)

                # This line printed to stdout in the original version as 
                # ofs was not openned as a file (i.e. no ofs.open(filename) line)
                oechem.OEWriteMolecule(ofs, mol)
            else:
                # Return an error if any problem with setting to pH 7.4
                sys.exit('ERROR: OESetNeutralpHModel returned False for ligand: {}'.format(ligand_prefix))
    
        ifs.close()
        ofs.close()
        os.rename(ligand_sdf_fpath, 
                  ligand_sdf_fpath[:-4]+'_before_pH_adjustment.sdf')
        os.rename(ligand_sdf_fpath[:-4]+'_pH74.sdf', ligand_sdf_fpath)

    #sp.call('obabel ' + ligand_sdf_fpath + ' -O ' + ligand_prefix + '.svg', shell=True)
 
    return 0


def oe_conformer_generation2(ligand_prefix_in, ligand_prefix_out, tauto_sp23=False, torsion_drive=True, box_cen=None, save_mol2=False, save_conf_isomer_ids=True):
    """
    Generate tautomers, enantiomers and conformers of the ligand. Ensure ligand 
    sdf file is in directory. 

    Parameters
    ----------
    ligand_prefix_in : str
        ligand file input name without file extension
    ligand_prefix_out : str
        ligand file output name without file extension
    tauto_sp23 : bool
        allow sp2 <-> sp3 conversion during tautomer enumeration
    torsion_drive : bool
        generate conformers with different torsion angles
    box_cen : list/array
        centre conformers at the point

    Returns
    -------
    filtered_by_pains : bool
        If True, the ligand was filtered out because it contains a PAINS fragment.
        If False, the ligand was not filtered out by the PAINS filter.

    """

    ligand_in_sdf_fpath = '{}.sdf'.format(ligand_prefix_in)
    ligand_out_sdf_fpath = '{}_confs.sdf'.format(ligand_prefix_out)
    #ligand_oeb_fpath = '{}.oeb'.format(ligand_prefix)
    
    # conformer options
    omega_opts = oeomega.OEOmegaOptions()
    omega_opts.SetEnergyWindow(20)
    omega_opts.SetMaxSearchTime(600)
    omega_opts.SetSearchForceField(7)
    omega_opts.SetRMSThreshold(0.5) 
    omega_opts.SetMaxConfs(1000) 

    # Turn off TorsionDrive to prevent search over torsion angles as this will be done in GNINA.  Still need to enumerate over stereoisomers and tautomers though.
    print('Use TorsionDrive: {}'.format(torsion_drive))
    omega_opts.SetTorsionDrive(torsion_drive)
    # Possible alternative way to turn off TorsionDrive:
    #torOpts = oeomega.OETorDriveOptions()
    #torOpts.ExceedsMaxRotors(0)
    #omega_opts.SetTorDriveOptions(torOpts)

    # Don't think the following option has any effect here since any undefinied stereochemistry will be defined by flipper which enumerates over undefined stereo:
    omega_opts.SetStrictStereo(True) # change to True is stereocenters are known
    omega = oeomega.OEOmega(omega_opts)

    # input molecule file
    ifs = oechem.oemolistream()
    ifs.open(ligand_in_sdf_fpath)

    # output molecule file
    #ofs = oechem.oemolostream()
    #ofs.open(ligand_oeb_fpath)
    ofs_sdf = oechem.oemolostream()
    ofs_sdf.SetFormat(oechem.OEFormat_SDF)
    ofs_sdf.open(ligand_out_sdf_fpath)

    # Output conformers in a mol2 file for showing in VMD
    if save_mol2:
        ofs_mol2 = oechem.oemolostream()
        ofs_mol2.SetFormat(oechem.OEFormat_MOL2)
        ofs_mol2.open(ligand_out_sdf_fpath[:-4]+'.mol2')
    
    if save_conf_isomer_ids:
        conf_isomer_ids = open(ligand_out_sdf_fpath[:-4]+'_conf_isomers.dat', 'w')
        conf_isomer_ids.write('conf_n,tauto,enant\n')

    if box_cen is not None:
        print('Translate to docking box centre: True')

    for mol in ifs.GetOEMols():
        
        # this function generates all conformers
        omega.Build(mol)

        ligand = oechem.OEMol(mol)

        # Optionally translate ligand to centre of docking box:
        if box_cen is not None:
            # Centre ligands at (0, 0, 0):
            oechem.OECenter(ligand)
            oe_box_cen = oechem.OEDoubleArray(3)
            for i in range(3):
                oe_box_cen[i] = box_cen[i]
            # Move ligands to docking box centre:
            oechem.OETranslate(ligand, oe_box_cen)

        #oechem.OEWriteMolecule(ofs, ligand)
        oechem.OEWriteMolecule(ofs_sdf, ligand)
        if save_mol2:
            oechem.OEWriteMolecule(ofs_mol2, ligand)

    ifs.close()
    #ofs.close()
    ofs_sdf.close()


def oe_gen_confs(ligand_smiles, ligand_prefix, tauto_sp23=False, ph=True, torsion_drive=True, box_cen=None):
    """
    Prepare the ligand and generate tautomers/enantiomers/conformers.

    Parameters
    ----------
    ligand_mol_fpath : str
        ligand file path (.mol or .sdf).

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # start_t = time.time()
    
    # ligand_prefix = ligand_mol_fpath[:-4]
    molecule_prep(ligand_smiles, ligand_prefix, ph)
    n_tauto, n_enant, n_disallowed, n_confs = \
    oe_conformer_generation(ligand_prefix, tauto_sp23=tauto_sp23, torsion_drive=torsion_drive, box_cen=box_cen)

    # Maybe keep all files just in case?
    #os.remove(ligand_oeb_fpath)


if __name__ == '__main__':
    parser = argparse.ArgumentParser() #exit_on_error=False)
    parser.add_argument('smi')
    parser.add_argument('molname')
    parser.add_argument('--tauto_sp23', action='store_true')
    parser.add_argument('--No_pH', action='store_true')
    parser.add_argument('--No_TorDrive', action='store_true')
    parser.add_argument('--center_x', type=float, nargs='?', default=None)
    parser.add_argument('--center_y', type=float, nargs='?', default=None)
    parser.add_argument('--center_z', type=float, nargs='?', default=None)
    # Convert into dictionary:
    args = vars(parser.parse_args())
    ligand_smiles = args.pop('smi')
    ligand_name = args.pop('molname')
    # Process some default options:
    if args.get('No_pH'):
        args['ph'] = False
    del args['No_pH']
    if args.get('No_TorDrive'):
        args['torsion_drive'] = False
    del args['No_TorDrive']
    if args['center_x'] and args['center_y'] and args['center_z']:
        args['box_cen']=[args.pop('center_x'), args.pop('center_y'), args.pop('center_z')]
    else:
        del args['center_x'], args['center_y'], args['center_z']
        args['box_cen']=None

    oe_gen_confs(ligand_smiles, ligand_name, **args) #, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
    #            ligand SMILES, ligand name, tauto_sp23, ph, torsion_drive, box_cen
    #oe_gen_confs(sys.argv[1], sys.argv[2]) #, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])