####################################################
options = {
    "level" : "R2SCAN-3C",
    "basis_set" : "",
    "solvent" : "phcf3",
    "step" : 0.02,

    "procs" : 16,
    "mem" : 4, # Memory per core, in GB
    "charge" : 0,
    "maxstep" : 0.1, # in Bohr atomic units (1au = 0.529177 A)
    }

epsilon = {
     "ch2cl2" : 9.04,
     "phcf3" : 9.18,
}
####################################################

import os 
import sys
from cclib.io import ccread
from numpy.linalg import norm

if len(sys.argv) == 1:
    print(f"\n  Makes an ORCA input file for a distance scan. Syntax:\n\n" +
           "  python mkorcascan.py filename.xyz i1 i2 target\n\n" + 
           "  filename.xyz: base name of input geometry file\n" +
           "  i1/i2: indices of atoms to scan\n" +
           "  target: desired distance between indices at end of scan.\n"
           )
    quit()


for option, value in options.items():
    print(f"--> {option} = {value}")

os.chdir(os.getcwd())
names = []

assert len(sys.argv) == 5, "Specify one input geometry file, two indices and one distance"

_, filename, i1, i2, d = sys.argv
i1, i2, d = int(i1), int(i2), float(d)
rootname = filename.split('.')[0]
coords = ccread(filename).atomcoords[0]
start_d = round(norm(coords[i1]-coords[i2]),2)
n_steps = round((d-start_d)/options["step"])
new_target_d = round(start_d + n_steps * options["step"], 2)
n_steps = abs(n_steps)+1

with open(f'{rootname}.inp', 'w') as f:
    s = f'''! {options["level"]} {options["basis_set"]} CPCM Opt

%pal
nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

%geom Scan
        B {i1} {i2} = {start_d}, {new_target_d}, {n_steps}
        end
      end

%cpcm
epsilon {epsilon[options["solvent"]]}
end

* xyzfile {options["charge"]} 1 {rootname}.xyz

'''

    f.write(s)
print(f'Written orca input file - {n_steps} ({options["step"]} A) steps from {start_d} A to {new_target_d} A')

# Convert all text files to Linux format
os.system(f'dos2unix {rootname}.inp')
os.system(f'dos2unix {rootname}.xyz')

print('Converted text files to Linux format.')