import argparse
import os
import textwrap

from run import run


def antegbsa(complex):
    """
    Run ante-MMPBSA.py to prepare topolopy files for MMPBSA.py

    Parameters
    ----------
    complex: str
        complex name

    Returns
    -------
    None (runs ante-MMPBSA.py)
    """

    ligandmask = complex.upper()

    if not os.path.exists('gbsa'):
        os.makedirs('gbsa')

    os.chdir('gbsa')
    os.system(f'rm -f {complex}-complex.parm7 {complex}-receptor.parm7 {complex}-ligand.parm7')

    run(f'ante-MMPBSA.py -p ../{complex}.parm7 -c {complex}-complex.parm7 -r {complex}-receptor.parm7 -l {complex}-ligand.parm7 --strip-mask=:WAT:CL:NA --ligand-mask=:{ligandmask} --radii=mbondi2')

    os.chdir('../')


def gbsa(complex):
    """
    Run MMPBSA.py

    Parameters
    ----------
    complex: str
        complex name

    Returns
    -------
    None (runs MMPBSA.py)
    """

    os.chdir('gbsa')

    string = textwrap.dedent('''\
    &general 
      interval=1, keep_files=1, netcdf=1,
    /   
    &gb 
      igb=5, saltcon=0.100,
    /
    ''')

    with open('gbsa.in', 'w') as f:
        f.write(string)

    run(f'MMPBSA.py -O -i gbsa.in -o {complex}_gbsa.dat -cp {complex}-complex.parm7 -rp {complex}-receptor.parm7 -lp {complex}-ligand.parm7 -y ../{complex}_equi_cent_strip.nc')

    os.chdir('../')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Run gbsa.py')

    # Required arguments
    parser.add_argument('-i','--input', help='Name of complexes to run functions on',required=True)

    # Optional arguments
    parser.add_argument('--cd', help='Name of directory to change into to run commands',required=False)
    parser.add_argument('-f','--functions', nargs='+', help='Function names to run',required=False)
    
    args = parser.parse_args()

    if args.cd:
        os.chdir(args.cd)

    complex = args.input

    if args.functions is None:
        antegbsa(complex)
        gbsa(complex)

    else:
        if 'antegbsa' in args.functions:
            antegbsa(complex)
        if 'gbsa' in args.functions:
            gbsa(complex)

    if args.cd:
        os.chdir('../')
