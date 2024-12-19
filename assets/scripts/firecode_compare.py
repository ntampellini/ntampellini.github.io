from firecode.utils import read_xyz, write_xyz, graphize
from prune import cl_similarity_refining
import numpy as np
import sys
import os
# from subprocess import getoutput

# e_thr = 20
atomnos = None
structures = []
energies = []
filenames = []

home = os.getcwd()

for name in sys.argv[1:]:

    folder = os.path.dirname(name)
    filename = os.path.basename(name)

    if folder != '':
        os.chdir(folder)

    print(f'Reading {filename}...', end='\r')
    mol = read_xyz(filename)
    filenames.append(filename)
    atomnos = atomnos if atomnos is not None else mol.atomnos
    for coords in mol.atomcoords:
        structures.append(coords)

    if folder != '':
        os.chdir(home)
    
    # energies_mol = getoutput(f'grep converged {name}').split('\n')
    # for line in energies_mol:
    #     energies.append(float(line.split()[0])*627.509608030593)

# energies = np.array(energies)
structures = np.array(structures, dtype=float)
before = len(structures)
print(f'--> XYZ similarity pruner: found {len(structures)} structures in {len(filenames)} files.')

###############

graph = graphize(structures[0], atomnos)

print('Removing similar structures...')
coords, _ = cl_similarity_refining(structures, atomnos, graph)
# energies = payload[0]
# structures, energies = zip(*sorted(zip(structures, energies), key=lambda x: x[1]))
# rel_energies = energies - np.min(energies)


outname = filenames[0].split('.')[0]+'.firecode.pooled.xyz'
with open(outname, "w") as f:
    # for i, (coord, energy) in enumerate(zip(structures, rel_energies)):
    for coord in structures:
        # if energy < e_thr:
        # title = f"Rel. E. = {energy:.3f} kcal/mol"
        write_xyz(coord, atomnos, f)

print(f"Wrote {len(structures)}/{before} structures to {outname}.")