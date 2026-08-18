from firecode.utils import read_xyz
import sys
from prism_pruner.algebra import get_alignment_matrix
import numpy as np

align_bond = (19, 36)
align_to = (0, 0, 1)

mol = read_xyz(sys.argv[1])
coords = mol.coords[0]
coords -= np.mean(coords, axis=0)

mat = get_alignment_matrix(np.array([align_to]), [coords[align_bond[1]]-coords[align_bond[0]]])
aligned_coords = (mat @ coords.T).T

mol.coords[0] = aligned_coords
outname = sys.argv[1].rstrip(".xyz") + "_aligned.xyz"
mol.to_xyz(outname)
print(f"Aligned {align_bond} to {align_to}.")