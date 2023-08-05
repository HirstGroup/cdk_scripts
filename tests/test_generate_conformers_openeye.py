import filecmp
import sys

sys.path.append('../')
from generate_conformers_openeye import *

def test2():

	oe_conformer_generation2('input/1-1', 'output/1-1', tauto_sp23=False, torsion_drive=False, box_cen=None, save_mol2=True, save_conf_isomer_ids=True)

def test_cycle():

	oe_conformer_generation2('input/cycle', 'output/cycle', tauto_sp23=False, torsion_drive=False, box_cen=None, save_mol2=True, save_conf_isomer_ids=True)

test_cycle()
