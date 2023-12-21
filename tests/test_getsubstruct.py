import filecmp
import sys
import textwrap

sys.path.append('../')
from util.getsubstruct import *

def test_get_substruct():

	atom_list = get_substruct('input/a33.sdf', '[#6]=[#6]-[#6]=O')

	assert atom_list == (0,1,2,3)


def test_has_pattern():

	assert has_pattern('input/a33.sdf', '[#6]=[#6]-[#6]=O') is True

def test_all():

	d = {'[#6]=[#6]-[#6]=O': (0,1,2,3), 'Cl[#6]-[#6]=O': (3,2,1,0), 'O=[#6]C#C': (3,2,1,0)}

	a = ['a33', 'a34', 'a35', 'a36', 'a44', 'a45', 'a46', 'a47', 'a48', 'a49']

	for lig in a:

		print(lig)

		for key, val in d.items():

			if has_pattern(f'input/{lig}.sdf', key):

				atom_list = get_substruct(f'input/{lig}.sdf', key)

				print(atom_list)

				assert atom_list == val

test_all()
