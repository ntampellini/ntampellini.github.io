from firecode.utils import read_xyz, write_xyz
import os
import sys
import argparse
from rich.traceback import install
from InquirerPy import inquirer
install(show_locals=True)

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="Input filename, in .xyz format.", action='store', default=None)
parser.add_argument("-o", "-O", "--outname", help="Output basename.", action='store', default=None)
parser.add_argument("-f", "--first", help="First structure index (1-based).", action='store', required=False, default=1)
parser.add_argument("-l", "--last", help="Last structure index (1-based).", action='store', required=False, default=-1)

args = parser.parse_args()
mol = read_xyz(args.inputfile)

if args.last == -1:
    args.last = len(mol.coords)

if inquirer.confirm(
    message=f"Split into {int(args.last)-int(args.first)+1} files? First name will be {args.outname}{args.first}.xyz",
    default=True,
    ):

    for c, coords in enumerate(mol.coords[int(args.first)-1:int(args.last)], start=int(args.first)):
        with open(f'{args.outname}{c}.xyz', 'w') as f:
            write_xyz(mol.atoms, coords, f)

    print(f'Wrote {int(args.last)-int(args.first)+1} files. The first is {args.outname}{args.first}.xyz, the last is {args.outname}{args.last}.xyz')

else:
    with open(f'{args.outname}.xyz', 'w') as f:
        for c, coords in enumerate(mol.coords[int(args.first)-1:int(args.last)], start=int(args.first)):
            write_xyz(mol.atoms, coords, f)

    print(f'Wrote {int(args.last)-int(args.first)+1} structures in {args.outname}.xyz')
