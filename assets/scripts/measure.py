import os 
import plotext as plt
import sys
from utils import read_xyz, norm_of, dihedral

plt.theme("pro")
plt.plotsize(100, 25)

def print_metadata(data):
    print()
    print(f"Avg {tag} = {sum(data)/len(data):8.3f} Å")
    print(f"Max {tag} = {max(data):8.3f} Å")
    print(f"Min {tag} = {min(data):8.3f} Å")
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
assert len(indices) in (2, 4)

cum_data = []

for filename in filenames:
    basename = filename.split(".")[0]   
    mol = read_xyz(basename+".xyz")

    if len(indices) == 2:
        tag = "distance (A)"
        i1, i2 = indices
        data = [norm_of(coords[i1]-coords[i2]) for coords in mol.atomcoords]
    else:
        tag = "dihedral angle (degrees)"
        i1, i2, i3, i4 = indices
        data = [dihedral(coords[indices]) for coords in mol.atomcoords]

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
