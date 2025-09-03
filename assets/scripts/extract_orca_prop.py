import os
import sys
from firecode.utils import read_xyz, write_xyz
from InquirerPy import inquirer

EH_TO_KCAL = 627.509608030593
BOHRS_TO_ANGSTROEM = 0.529177

def extract_geom_from_propfile(filename):

    xyzname = filename[:-4] + '_extracted.xyz'
    if xyzname in os.listdir():
        print(f'--> Skipping {xyzname}, already extracted.')
    else:

        with open(filename, 'r') as f:

            s = ''
            n_atoms = 0

            while True:
                line = f.readline()
                if '&CartesianCoordinates' in line:

                    while True:
                        line = f.readline()

                        if "$End" in line:
                            break

                        s += line
                        n_atoms += 1

                if '&FINALEN' in line:
                    energy = float(line.split()[-6])
                    break

                if not line:
                    break

        with open(xyzname, 'w') as f:
            
            f.write(str(n_atoms))
            f.write(f'\nEnergy (Eh) = {energy}\n')            
            f.write(s)
            print(f'--> Extracted geom to {xyzname}')

    mol = read_xyz(xyzname)
    coords = mol.atomcoords[0] * BOHRS_TO_ANGSTROEM
    atomnos = mol.atomnos

    return coords, atomnos

if __name__ == '__main__':

    # tgt = sys.argv[1]
    # os.chdir(tgt)

    structures = []
    # for filename in os.listdir():
    for filename in sys.argv[1:]:
        if filename.endswith('.property.txt'):
            try:
                coords, atomnos = extract_geom_from_propfile(filename)
                structures.append(coords)
            except Exception:
                pass

    outname = inquirer.text(
        message='Cumulative output filename?',
        filter=lambda x: x if x.endswith(".xyz") else x + '.xyz',
    ).execute()

    with open(outname, 'w') as f:
        for structure in structures:
            write_xyz(structure, atomnos, f)

    print(f'--> Printed {len(structures)} structures to {outname}')

        
