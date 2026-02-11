# Set of utility functions for the scripts at https://ntampellini.github.io/scripts/

import re

from firecode.pt import pt
from firecode.utils import read_xyz
from firecode.units import EH_TO_KCAL
from subprocess import getoutput

def d_min_bond(e1, e2, factor=1.2):
    return factor * (pt[e1].covalent_radius + pt[e2].covalent_radius)

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
                electrons += getattr(pt, element).number
            except AttributeError:
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
    a1, a2 = pt[mol.atomnos[i1]], pt[mol.atomnos[i2]]
    cr1 = a1.covalent_radius
    cr2 = a2.covalent_radius

    est_d = round(factor * (cr1 + cr2), 2)

    if verbose:
        print(f'--> Estimated TS d({a1}-{a2}) = {est_d} Å')
        
    return est_d


def read_xyz_energies(filename, verbose=True):
    '''
    Read energies from a .xyz file. Returns None or an array of floats (in Hartrees).
    '''
    energies = None

    # get lines right after the number of atom, which should contain the energy
    comment_lines = getoutput(f'grep -A1 "^[[:space:]]*[0-9]\\+$" {filename} | grep -v "^[[:space:]]*[0-9]\\+$" | grep -v "^--$"').split("\n")

    if len(comment_lines[0].split()) == 1:
        if set(comment_lines[0].split()[0]).issubset('0123456789.-'):
            # only one energy found with no UOM, assume it's in Eh
            energies = [float(e.split()[0].strip()) for e in comment_lines]

            if verbose:
                print(f'--> Read {len(energies)} energies from {filename} (single number, no UOM: assuming Eh units).')

        else:
            if verbose:
                print(f'--> Could not parse energies for {filename} - skipping.')

    else:
        # multiple energies found, parse units
        hartree_matches = re.findall(r'-*\d+\.\d+\sEH', comment_lines[0].upper())
        kcal_matches = re.findall(r'-*\d+\.\d+\sKCAL/MOL', comment_lines[0].upper())
        number_matches = re.findall(r'-*\d+\.\d+', comment_lines[0])

        if hartree_matches:
            energies = [float(re.findall(r'-*\d+\.\d+\sEH', e.upper())[0].split()[0].strip()) for e in comment_lines]
            if verbose:
                print(f'--> Read {len(comment_lines)} energies from {filename} (first number followed by Eh units).')

        elif kcal_matches:
            energies = [float(re.findall(r'-*\d+\.\d+\sKCAL/MOL', e.upper())[0].split()[0].strip())/EH_TO_KCAL for e in comment_lines]
            if verbose:
                print(f'--> Read {len(comment_lines)} energies from {filename} (first number followed by kcal/mol units).')
    
        # last resort, parse the first thing that looks like an energy and assume it's in Eh
        elif number_matches:
            energies = [float(re.findall(r'-*\d+\.\d+', e)[0].strip()) for e in comment_lines]
            if verbose:
                print(f'--> Read {len(comment_lines)} energies from {filename} (first number, no UOM: assuming Eh units).')

        else:
            if verbose:
                print(f'--> Could not parse energies for {filename} - skipping.')

    return energies
