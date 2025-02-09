import sys
import os

import numpy as np
from firecode.algebra import vec_angle
from firecode.utils import read_xyz, write_xyz

from InquirerPy import inquirer

# setting ids for the key atoms

indices_dict = {
    "bzh_arm" : (107, 106, 108, 119),
    "anth" : (104, 103, 105, 116),
    "tBuNHCO" : (97, 96, 98, 109),
    "dPhg" : (102, 101, 103, 114),
    "homoBzh" : (102, 103, 114, 101),
}

set_name = inquirer.select(
        message="Select which stereochemistry indices set to use:",
        choices=list(indices_dict.keys()),
        default="bzh_arm",
    ).execute()

key_atom, sub1, sub2, sub3 = indices_dict[set_name]
thr = 45
def_chir = 'S'

opposite_chirality = {'R':'S', 'S':'R'}

def get_absolute(coords, index=key_atom, sub1=sub1, sub2=sub2, sub3=sub3, thr=thr, def_chir=def_chir):
    '''
    '''
    v1 = coords[sub1] - coords[index]
    v2 = coords[sub2] - coords[index]
    v3 = coords[sub3] - coords[index]
    v_cross = np.cross(v1, v2)
    angle = vec_angle(v3, v_cross)

    if angle < thr:
        chirality = def_chir
    elif np.abs(angle-180) < thr:
        chirality = opposite_chirality[def_chir]
    else:
        chirality = f'Unknown - d{def_chir}={angle:5.2f}°, d{opposite_chirality[def_chir]}={np.abs(angle-180):5.2f}°'

    return chirality

if __name__ == '__main__':    

    if "thr" in [kw.split("=")[0] for kw in sys.argv]:
        thr = next((kw.split("=")[-1] for kw in sys.argv if "thr=" in kw))
        sys.argv.remove(f"thr={thr}")

    print('\n--> Absolute configuration reader script')
    print(f'  key: {key_atom}')
    print(f'  s1: {sub1}')
    print(f'  s2: {sub2}')
    print(f'  s3: {sub3}')
    print(f'  thr: {thr}')
    print(f'  default chirality: {def_chir}\n')

    extract_abs = None
    if "R" in sys.argv:
        sys.argv.remove("R")
        extract_abs = 'R'
    if "S" in sys.argv:
        sys.argv.remove("S")
        extract_abs = 'S'

    if "x" in [kw.split("=")[0] for kw in sys.argv]:
        outname = next((kw.split("=")[-1] for kw in sys.argv if "=" in kw))
        try:
            sys.argv.remove(f"x={outname}")
        except ValueError:
            raise ValueError(f'Trying to delete non-existent\"x={outname}\"')
    else:
        outname = None

    coords_to_extract = []
    total = 0
    for filename in sys.argv[1:]:
        mol = read_xyz(filename)
        for c, coords in enumerate(mol.atomcoords):
            total += 1
            chirality = get_absolute(coords)
            print(f'{filename:20s}, conf {c:3} - {chirality}')

            if extract_abs is not None:
                if chirality == extract_abs:
                    coords_to_extract.append((filename, coords, mol.atomnos))

    if extract_abs is not None:

        print(f'--> {len(coords_to_extract)}/{total} structures have the desired {extract_abs} configuration.')

        if outname is not None:

            if os.path.basename(outname) == '':
                print(f'--> Copying {len(coords_to_extract)} structures to same-name new files in {os.path.dirname(outname)}')

                for filename, coords, atomnos in coords_to_extract:

                    outname = os.path.join(os.path.dirname(outname), filename)[:-4] + '.xyz'
                    with open(outname, "w") as f:
                        title = f"{filename} - {extract_abs} configuration at atom {key_atom}"
                        write_xyz(coords, atomnos, f, title=title)

            else:
                with open(outname, "w") as f:
                    for filename, coords, atomnos in coords_to_extract:
                        title = f"{filename} - {extract_abs} configuration at atom {key_atom}"
                        write_xyz(coords, atomnos, f, title=title)

                print(f"Wrote {len(coords_to_extract)} structures to {outname}.")

