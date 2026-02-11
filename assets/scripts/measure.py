import os 
import numpy as np
import plotext as plt
import sys
from prism_pruner.algebra import dihedral
from firecode.algebra import point_angle
from firecode.utils import read_xyz

plt.theme("pro")
plt.plotsize(100, 25)

def print_metadata(data):
    s = '-'.join([str(mol.atoms[i])+'('+str(i)+')' for i in indices])
    print(f'\n--> Showing {s} {tag}.')
    print(f"  Min = {min(data):8.3f} Å")
    print(f"  Max = {max(data):8.3f} Å")
    print(f"  Avg = {sum(data)/len(data):8.3f} Å")
    print()

if len(less_than:=[kw for kw in sys.argv if "<" in kw]) > 0:
    sys.argv.remove(less_than[0])
    less_than = less_than[0][1:]
    print(f'--> Only showing values smaller than {less_than}')

if len(more_than:=[kw for kw in sys.argv if ">" in kw]) > 0:
    sys.argv.remove(more_than[0])
    more_than = more_than[0][1:]
    print(f'--> Only showing values greater than {less_than}')

basenames = []
indices = [arg for arg in sys.argv if "." not in arg]
filenames = set(sys.argv[1:]) - set(indices)

# determine style of results
indices = [int(i) for i in indices]
assert len(indices) in (2, 3, 4)

cum_data = []

for filename in filenames:
    basename = filename[:-4]
    mol = read_xyz(basename+".xyz")

    if len(indices) == 2:
        tag = "distance (A)"
        i1, i2 = indices
        data = [np.linalg.norm(coords[i1]-coords[i2]) for coords in mol.coords]
    
    elif len(indices) == 3:
        tag = "planar angle (degrees)"
        i1, i2, i3 = indices
        data = [point_angle(*coords[indices]) for coords in mol.coords]

    else:
        tag = "dihedral angle (degrees)"
        i1, i2, i3, i4 = indices
        data = [dihedral(coords[indices]) for coords in mol.coords]

    if len(filenames) == 1:
        plt.simple_bar(data, width=50)
        plt.xlabel("Structure #")
        plt.ylabel(f"{indices} {tag}")
        plt.show()
        print_metadata(data)
        sys.exit()

    else:
        if not less_than and not more_than:
            cum_data.append(data[-1])
            basenames.append(basename)

        elif less_than and data[-1] < float(less_than):
            cum_data.append(data[-1])
            basenames.append(basename)

        elif more_than and data[-1] > float(more_than):
            cum_data.append(data[-1])
            basenames.append(basename)

plt.simple_bar(basenames, cum_data, width=50)
plt.xlabel("Files")
plt.ylabel(f"{indices} {tag}")
print_metadata(cum_data)

plt.show()
