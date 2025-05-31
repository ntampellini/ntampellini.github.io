# Set of utility functions for the scripts at https://ntampellini.github.io/scripts/

import os

import numpy as np
from cclib.io import ccread
from networkx import from_numpy_matrix, set_node_attributes
from periodictable import core, covalent_radius, mass

for pt_n in range(5):
    try:
        pt = core.PeriodicTable(table=f"H={pt_n+1}")
        covalent_radius.init(pt)
        mass.init(pt)
    except ValueError:
        continue
    break

def d_min_bond(e1, e2, factor=1.2):
    return factor * (pt[e1].covalent_radius + pt[e2].covalent_radius)

class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)

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

def point_angle(p1, p2, p3):
    '''
    Returns the planar angle between three points in space, in degrees.
    '''
    return np.arccos(np.clip(norm(p1 - p2) @ norm(p3 - p2), -1.0, 1.0))*180/np.pi

def norm_of(vec):
    '''
    Returns the norm of the vector.
    Faster than Numpy version, but 
    only compatible with 3D vectors.
    '''
    return ((vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))**0.5

def vec_angle(v1, v2):
    v1_u = norm(v1)
    v2_u = norm(v2)
    return np.arccos(clip(np.dot(v1_u, v2_u), -1.0, 1.0))*180/np.pi

def clip(n, lower, higher):
    '''
    jittable version of np.clip for single values
    '''
    if n > higher:
        return higher
    elif n < lower:
        return lower
    else:
        return n

def point_angle(p1, p2, p3):
    return np.arccos(np.clip(norm(p1 - p2) @ norm(p3 - p2), -1.0, 1.0))*180/np.pi

def norm(vec):
    '''
    Returns the normalized vector.
    Reasonably faster than Numpy version.
    Only for 3D vectors.
    '''
    return vec / np.sqrt((vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]))

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

def all_dists(arr1, arr2):
    '''
    Returns an array with all the distances between vectors
    in the two arrays.

    '''
    assert arr1.shape == arr2.shape
    output = np.zeros((arr1.shape[0], arr1.shape[0]), dtype=float)

    for i1, v1 in enumerate(arr1):
        for i2, v2 in enumerate(arr2):
            if i2 > i1:
                output[i1][i2] = norm_of(v1-v2)

    return output

def graphize(coords, atomnos, mask=None):
    '''
    :params coords: atomic coordinates as 3D vectors
    :params atomnos: atomic numbers as a list
    :params mask: bool array, with False for atoms
    to be excluded in the bond evaluation
    :return connectivity graph
    
    '''

    mask = np.array([True for _ in atomnos], dtype=bool) if mask is None else mask

    matrix = np.zeros((len(coords),len(coords)))
    for i, _ in enumerate(coords):
        for j in range(i,len(coords)):
            if mask[i] and mask[j]:
                if norm_of(coords[i]-coords[j]) < d_min_bond(atomnos[i], atomnos[j]):
                    matrix[i][j] = 1

    graph = from_numpy_matrix(matrix)
    set_node_attributes(graph, dict(enumerate(atomnos)), 'atomnos')

    return graph
