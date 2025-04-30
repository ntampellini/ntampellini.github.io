import os as op_sys
import sys
from subprocess import getoutput

maxlen = max(len(f) for f in sys.argv[1:])

print(f"{'FILENAME':{maxlen}s} {'#ATOMS':3s} {'CHARGE':2s}")
print(f"--------------------------------------------------------")
for filename in sys.argv[1:]:
    if not "_trj" in filename and not "Compound" in filename:
        atoms = getoutput(f"head {filename} -n 1")
        charge = str(filename.count("+") - filename.count("-"))
        print(f"{filename:{maxlen}s} {atoms:3s} {charge:>2s}")