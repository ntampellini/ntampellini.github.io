from firecode.utils import read_xyz, write_xyz
import os
import sys
import argparse
from rich.traceback import install
install(show_locals=True)

# if len(sys.argv) == 1:
#     print('usage: split.py source rootname')

parser = argparse.ArgumentParser()
parser.add_argument("inputfile", help="Input filename, in .xyz format.", action='store', default=None)
parser.add_argument("out_basename", help="Output basename.", action='store', default=None)
parser.add_argument("-f", "--first", help="First structure index (1-based).", action='store', required=False, default=1)
parser.add_argument("-l", "--last", help="Last structure index (1-based).", action='store', required=False, default=-1)

args = parser.parse_args()
mol = read_xyz(args.inputfile)

for c, coords in enumerate(mol.atomcoords[int(args.first)-1:int(args.last)], start=1):
    with open(f'{args.out_basename}{c}.xyz', 'w') as f:
        write_xyz(coords, mol.atomnos, f)

print(f'Wrote {int(args.last)-int(args.first)+1} files. The first is {args.out_basename}1.xyz')
