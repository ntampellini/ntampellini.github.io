import os 
import numpy as np
import plotext as plt
import sys
from prism_pruner.algebra import dihedral
from firecode.algebra import point_angle
from firecode.utils import read_xyz
from prism_pruner.graph_manipulations import graphize

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
    less_than = float(less_than[0][1:])
    print(f'--> Only showing values smaller than {less_than}')
else:
    less_than = None

if len(more_than:=[kw for kw in sys.argv if ">" in kw]) > 0:
    sys.argv.remove(more_than[0])
    more_than = float(more_than[0][1:])
    print(f'--> Only showing values greater than {less_than}')
else:
    more_than = None

def check_value(value: float) -> bool:
    if less_than is not None:
        if value > less_than:
            return False
        
    if more_than is not None:
        if value < more_than:
            return False
        
    return True

basenames = []
indices = [arg for arg in sys.argv if "." not in arg]
filenames = list(set(sys.argv[1:]) - set(indices))

if all(i.isdigit() for i in indices):

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
            if check_value(data[-1]):
                cum_data.append(data[-1])
                basenames.append(basename)

    plt.simple_bar(basenames, cum_data, width=50)
    plt.xlabel("Files")
    plt.ylabel(f"{indices} {tag}")
    print_metadata(cum_data)
    plt.show()

elif all(i.isalpha() for i in indices):

    if len(indices) != 2:
        raise NotImplementedError(f'Please only specify two symbols. len({indices}) == {len(indices)}')

    if len(filenames) > 1:
        print(f'--> More than one filename specified: will only work on the first.')

    mol = read_xyz(filenames[0])
    symbols = (indices[0], indices[1])

    graph = graphize(mol.atoms, mol.coords[0])
    bond_symbols = {(i, j): (mol.atoms[i], mol.atoms[j]) for i, j in graph.edges}

    matches = dict()
    for (i, j), (e1, e2) in bond_symbols.items():
        if (e1, e2) == symbols or (e2, e1) == symbols:
            d = np.linalg.norm(mol.coords[0][i]-mol.coords[0][j])
            if check_value(d):
                    matches[(i, j)] = (d, 'bonded')

    for i1, e1 in enumerate(mol.atoms):
        for i2, e2 in enumerate(mol.atoms):
            if i2 > i1 and (e1, e2) == symbols or (e2, e1) == symbols:
                if (i1, i2) not in matches.keys() and (i2, i1) not in matches.keys():
                    d = np.linalg.norm(mol.coords[0][i1]-mol.coords[0][i2])
                    if check_value(d):
                        matches[(i1, i2)] = (d, 'not bonded')

    if matches:

        # sort by distance
        matches = dict(sorted(matches.items(), key=lambda x: x[1][0]))

        print()

        for (i, j), (d, bonded) in matches.items():
            print(f'  {mol.atoms[i]:2}({i:>3})-{mol.atoms[j]:2}({j:>3}) : {d:.2f} Å ({bonded})')
        
        print()

    else:
        print(f'No matches found for {indices[0]}-{indices[1]} distances.')

else:
    print('Input not understood. Only specify .xyz filename(s) and indices or chemical symbols.')
