from firecode.utils import read_xyz, write_xyz, graphize
from prune import cl_similarity_refining
import numpy as np
import os
from subprocess import getoutput

e_thr = 10
atomnos = None
structures = []
energies = []
filenames = []
for name in os.listdir():
    if 'finalensemble.xyz' in name:
        filenames.append(name)
        print(f'Reading {name}...', end='\r')
        mol = read_xyz(name)
        atomnos = atomnos if atomnos is not None else mol.atomnos
        for coords in mol.atomcoords:
            structures.append(coords)
        
        energies_mol = getoutput(f'grep converged {name}').split('\n')
        for line in energies_mol:
            energies.append(float(line.split()[0])*627.509608030593)

energies = np.array(energies)
structures = np.array(structures, dtype=float)
before = len(structures)
print(f'--> GOAT extractor: found {len(structures)} structures in {len(filenames)} files.')

###############

graph = graphize(structures[0], atomnos)

print('Removing similar structures...')
coords, payload = cl_similarity_refining(structures, atomnos, graph, payload=[energies])
energies = payload[0]
structures, energies = zip(*sorted(zip(structures, energies), key=lambda x: x[1]))
rel_energies = energies - np.min(energies)


outname = filenames[0].split('.')[0]+'.goat.pooled.xyz'
with open(outname, "w") as f:
    for i, (coord, energy) in enumerate(zip(structures, rel_energies)):
        if energy < e_thr:
            title = f"Rel. E. = {energy:.3f} kcal/mol"
            write_xyz(coord, atomnos, f, title=title)

print(f"Wrote {len(structures)}/{before} structures to {outname}.")