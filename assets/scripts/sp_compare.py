from firecode.utils import read_xyz, write_xyz
import sys
import os
from subprocess import getoutput

def read_gcorr(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        if "gcorr(Eh)" in line:
            return float(line.split()[2])
        
def read_sp_ee(outname):
    return float(getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {outname}').split()[4])

extract = False
if "x" in [kw.split("=")[0] for kw in sys.argv]:
    outname = next((kw.split("=")[-1] for kw in sys.argv if "x" in kw))
    sys.argv.remove(f"x={outname}")
    extract = True

all = False
if "all" in sys.argv:
    sys.argv.remove("all")
    all = True

def main(argvs, outname):

    data = dict()

    for filename in argvs[1:]:
        basename = filename.split(".")[0]

        gcorr = read_gcorr(basename+'.xyz')
        sp_ee = read_sp_ee(basename+'.out')
        
        g = sp_ee + gcorr

        data[basename] = {"G(Eh)":g}
        data[basename]['EE(Eh)'] = sp_ee
        data[basename]['gcorr(Eh)'] = gcorr

    min_e = min([entry["G(Eh)"] for entry in data.values()])

    for name, d in data.items():
        d["Rel. G (kcal/mol)"] = (d["G(Eh)"] - min_e) * 627.509608030592

    print("Name                Rel. G (kcal/mol)")
    print("-------------------------------------")
    for basename, d in sorted(data.items(), key=lambda item: item[1]["G(Eh)"]):
        print(f"{basename:20}{d['Rel. G (kcal/mol)']:.2f}")

    if extract:
        with open(outname, "w") as f:
            for basename, d in sorted(data.items(), key=lambda item: item[1]["G(Eh)"]):
                mol = read_xyz(basename+".xyz")
                write_xyz(mol.atomcoords[0], mol.atomnos, f, title=f"SP_EE(Eh) = {d['EE(Eh)']:.8f}, gcorr(Eh) = {d['gcorr(Eh)']:.8f}, G(Eh) = {d['G(Eh)']:.8f}, Rel. G(kcal/mol) = {d['Rel. G (kcal/mol)']:.2f}")

        print(f'Written {len(data.items())} structures to {outname}.')

if __name__ == '__main__':

    if all:
        extract = True
        for basename in {filename.split("_")[0] for filename in os.listdir() if filename.endswith('.xyz')}:
            print(f'Extracting {basename}')
            outname = basename + "_sp.xyz"

            files = getoutput(f'ls {basename}*[^sp].xyz').split()

            main([None, *files], outname)

    else:
        main(sys.argv, outname)