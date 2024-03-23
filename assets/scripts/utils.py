# Set of utility functions for the scripts at https://ntampellini.github.io/scripts/

import numpy as np
from periodictable import core, covalent_radius, mass
from cclib.io import ccread

for pt_n in range(5):
    try:
        pt = core.PeriodicTable(table=f"H={pt_n+1}")
        covalent_radius.init(pt)
        mass.init(pt)
    except ValueError:
        continue
    break

def cycle_to_dihedrals(cycle):
    '''
    '''
    dihedrals = []
    for i in range(len(cycle)):

        a = cycle[i % len(cycle)]
        b = cycle[(i+1) % len(cycle)]
        c = cycle[(i+2) % len(cycle)]
        d = cycle[(i+3) % len(cycle)]
        dihedrals.append([a, b, c, d])
    return dihedrals

def get_exocyclic_dihedrals(graph, cycle):
    '''
    '''
    exo_dihs = []
    for index in cycle:
        for exo_id in neighbors(graph, index):
            if exo_id not in cycle:
                dummy1 = next(i for i in cycle if i not in (exo_id, index) and i in neighbors(graph, index))
                dummy2 = next(i for i in cycle if i not in (exo_id, index, dummy1) and i in neighbors(graph, dummy1))
                exo_dihs.append([exo_id, index, dummy1, dummy2])

    return exo_dihs 

def neighbors(graph, index):
    # neighbors = list([(a, b) for a, b in graph.adjacency()][index][1].keys())
    neighbors = list(graph.neighbors(index))
    if index in neighbors:
        neighbors.remove(index)
    return neighbors

def dihedral(p):
    '''
    Returns dihedral angle in degrees from 4 3D vecs
    Praxeolitic formula: 1 sqrt, 1 cross product
    
    '''
    p0 = p[0]
    p1 = p[1]
    p2 = p[2]
    p3 = p[3]

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= norm_of(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    
    return np.degrees(np.arctan2(y, x))

def norm_of(vec):
    '''
    Returns the norm of the vector.
    Faster than Numpy version, but 
    only compatible with 3D vectors.
    '''
    return ((vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))**0.5

def write_xyz(coords:np.array, atomnos:np.array, output, title='temp'):
    '''
    output is of _io.TextIOWrapper type

    '''
    assert atomnos.shape[0] == coords.shape[0]
    assert coords.shape[1] == 3
    string = ''
    string += str(len(coords))
    string += f'\n{title}\n'
    for i, atom in enumerate(coords):
        string += '%s     % .6f % .6f % .6f\n' % (pt[atomnos[i]].symbol, atom[0], atom[1], atom[2])
    output.write(string)

def read_xyz(filename):
    '''
    Wrapper for ccread. Raises an error if unsuccessful.

    '''
    mol = ccread(filename)
    assert mol is not None, f'Reading molecule {filename} failed - check its integrity.'
    return mol