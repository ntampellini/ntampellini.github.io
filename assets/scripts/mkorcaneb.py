#%%

import sys
from tscode.utils import read_xyz, write_xyz
import numpy as np

assert len(sys.argv) == 2, "Usage: mkorcaneb.py mep_guess.xyz"

scan_geoms = sys.argv[1]
mol = read_xyz(sys.argv[1])
n_scan_geoms = len(mol.atomcoords)
n_images = 9
print(f"Scan file has {n_scan_geoms} geoms, extracting {n_images}:")

#%%
guess_ids = [int((i+1) * n_scan_geoms / n_images)-1 for i in np.arange(n_images)]
guess_ids[0] = 0
print(guess_ids)

with open("mep_guess.xyz", "w") as f:
    for guess_id in guess_ids:
        write_xyz(mol.atomcoords[guess_id], mol.atomnos, f)
#%%
with open("mep_guess.xyz", "r") as f:
    atoms_string = f.readline()
    s = atoms_string

    while True:
        line = f.readline()

        if not line:
            break

        if line == atoms_string:
            s += ">\n"

        s += line

with open("mep_guess.allxyz", "w") as f:
    f.write(s)

with open("end.xyz", "w") as f:
    write_xyz(mol.atomcoords[-1], mol.atomnos, f)

with open("start.xyz", "w") as f:
    write_xyz(mol.atomcoords[0], mol.atomnos, f)



inpname = "neb.inp"

s = f'''! R2SCAN-3c CPCM NEB-TS
! Defgrid3

%pal
  nprocs 18
end

%maxcore 5000

%geom
  MaxStep 0.05
end

%cpcm
  epsilon 4.81
end

%neb
  product "end.xyz"
  NImages 9
  Restart_ALLXYZFile "mep_guess.allxyz"
  Free_End true
  Stepsize 1 # multiplicative factor, not abs distance
  Maxiter 500
end

* xyzfile 0 1 start.xyz
'''

with open(inpname, "w") as f:
    f.write(s)