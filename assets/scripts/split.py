from firecode.utils import read_xyz, write_xyz
import os
import sys
import numpy as np
import argparse
from rich.traceback import install
from InquirerPy import inquirer
install(show_locals=True)

parser = argparse.ArgumentParser()
parser.add_argument("inputfiles", help="Input filename(s), in .xyz format.", action='store', nargs="*", default=None)
parser.add_argument("-o", "-O", "--outname", help="Output filename (\"mol.xyz\", single file) or basename (\"mol\", multiple files).", action='store', default=None)
parser.add_argument("-f", "--first", help="First structure index (1-based).", action='store', required=False, default=1)
parser.add_argument("-l", "--last", help="Last structure index (1-based).", action='store', required=False, default=-1)

args = parser.parse_args()

all_coords = []
atoms = None
for file in args.inputfiles:
    mol = read_xyz(file)

    if atoms is not None:
        assert np.all(mol.atoms == atoms), "All files must have the same atoms!"
        atoms = mol.atoms

    all_coords.extend(mol.coords)

print(f"--> Collected {len(all_coords)} total structures.")

if args.last == -1:
    args.last = len(all_coords)

if args.first == -1:
    args.first = len(all_coords)

if args.outname.endswith(".xyz"):

    with open(f'{args.outname}', 'w') as f:
        for c, coords in enumerate(all_coords[int(args.first)-1:int(args.last)], start=int(args.first)):
            write_xyz(mol.atoms, coords, f)

    print(f'Wrote {int(args.last)-int(args.first)+1} structures in {args.outname}')

else:
    for c, coords in enumerate(all_coords[int(args.first)-1:int(args.last)], start=int(args.first)):
        with open(f'{args.outname}{c}.xyz', 'w') as f:
            write_xyz(mol.atoms, coords, f)

    print(f'Wrote {int(args.last)-int(args.first)+1} files. The first is {args.outname}{args.first}.xyz, the last is {args.outname}{args.last}.xyz')
