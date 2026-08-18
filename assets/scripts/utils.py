# Set of utility functions for the scripts at https://ntampellini.github.io/scripts/

import re

from firecode.pt import pt
from firecode.utils import read_xyz
from firecode.units import EH_TO_KCAL
from subprocess import getoutput

USE_COLOR = True

class tcolors:
    HEADER = '\033[95m' if USE_COLOR else ""
    BLUE = '\033[94m' if USE_COLOR else ""
    CYAN = '\033[96m' if USE_COLOR else ""
    GREEN = '\033[38;2;65;165;165m' if USE_COLOR else ""
    YELLOW = '\033[93m' if USE_COLOR else ""
    RED = '\033[38;2;232;111;136m' if USE_COLOR else ""
    ENDC = '\033[0m' if USE_COLOR else ""
    BOLD = '\033[1m' if USE_COLOR else ""
    UNDERLINE = '\033[4m' if USE_COLOR else ""

def d_min_bond(e1, e2, factor=1.2):
    return factor * (pt.covalent_radius(e1) + pt.covalent_radius(e2))

def multiplicity_check(rootname, charge, multiplicity=1) -> bool:
    '''
    Returns True if the multiplicity and the nuber of
    electrons are one odd and one even, and vice versa.

    '''

    electrons = 0
    for line in getoutput(f'cat {rootname}.xyz').splitlines():
        parts = line.split()
        if len(parts) == 4:
            try:
                element = parts[0]
                electrons += pt.number(element)
            except KeyError:
                pass

    electrons -= charge
    
    return (multiplicity % 2) != (electrons % 2)

def get_ts_d_estimate(filename, indices, factor=1.35, verbose=True):
    '''
    Returns an estimate for the distance between two
    specific atoms in a transition state, by multipling
    the sum of covalent radii for a constant.
    
    '''
    mol = read_xyz(filename)
    i1, i2 = indices
    s1, s2 = mol.atoms[i1], mol.atoms[i2]
    cr1 = pt.covalent_radius(s1)
    cr2 = pt.covalent_radius(s2)

    est_d = round(factor * (cr1 + cr2), 2)

    if verbose:
        print(f'--> Estimated TS d({s1}-{s2}) = {est_d} Å')
        
    return est_d


