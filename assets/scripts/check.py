import plotext as plt
from cclib.io import ccread
import sys
from numpy.linalg import norm
import numpy as np
from subprocess import getoutput
import os
from utils import dihedral, write_xyz

plt.theme("pro")
plt.plotsize(100,25)
plot_distance = False
scan = False
neb = False
indices = None
zoom = False

if len(sys.argv) == 1:
    print("\n  Check running and completed ORCA jobs and prints/plots relevant information. " +
          "For scans, extracts local maxima structures. Syntax:\n\n" +
          "  python check.py filename[.out] [i1] [i2] [zoom]\n\n" + 
          "  filename: base name of input/output file (with or without extension)\n" +
          "  i1/i2: optional, plot i1-i2 distance during optimization steps\n" +
          "  zoom: optional, only shows plot of last 20 iterations.\n")
    quit()

# move to target directory and get basename
remote_dir = os.path.dirname(sys.argv[1])
if remote_dir != "":
    os.chdir(remote_dir)
    sys.argv[1] = os.path.basename(sys.argv[1])

if "zoom" in sys.argv:
    sys.argv.remove("zoom")
    zoom = True

if len(sys.argv) > 2:
    indices = [int(i) for i in sys.argv[2:]]
    sys.argv = sys.argv[:2]

_, rootname = sys.argv
rootname = rootname.split(".")[0]
inpname = f"{rootname}.inp"
files = os.listdir()

if inpname in files:
    with open(inpname, "r") as f:
        while True:

            line = f.readline()
            if not line:
                break

            if "SCAN" in line.upper():
                scan = True
                frags = f.readline().upper().split()

                if frags[0] == "B":
                    scantype = "distance"
                    scan_indices = tuple([int(i) for i in frags[1:3]])

                elif frags[0] == "D":
                    scantype = "dihedral"
                    scan_indices = tuple([int(i) for i in frags[1:5]])

                else:
                    scan = False

            if "NEB" in line.upper():
                neb = True

filename = rootname + "_trj.xyz"

if filename in files and indices is not None:

    mol = ccread(filename)

    if len(indices) == 2:
        y = [norm(c[indices[0]]-c[indices[1]]) for c in mol.atomcoords]
        tag = "distance (A)"

    elif len(indices) == 4:
        y = [dihedral(c[indices]) for c in mol.atomcoords]
        tag = "dihedral angle (Degrees)"

    else:
        raise Exception("Provide 2 or 4 indices")

    plot = plt.plot(y)
    plt.xlabel("Iteration #")
    plt.ylabel(f"{indices} {tag}")
    plt.show()
    print("\n"+"_"*100+"\n")

propname = f"{rootname}_property.txt"

if propname in files:
    
    if not neb:
        with open(propname, 'r') as f:
            energies = []
            while True:
                line = f.readline()
                if "Total DFT Energy" in line:
                    energies.append(float(line.split()[6]))

                if not line:
                    break

    else:
        energies = []
        lines = getoutput(f"grep \"Starting iterations:\" {rootname}.out -A 500 ").splitlines()
        # energies = [float(line.split()[3]) for line in lines[4:]]
        for line in lines[4:]:
            try:
                assert line.split()[0] != "Convergence"
                energies.append(float(line.split()[3]))

            except (IndexError, ValueError, AssertionError):
                continue
   
    if len(energies) > 0:
        last_E = energies[-1]
        energies = np.array(energies)
        energies -= np.min(energies) if not neb else 0
        energies *= 627.5096080305927
        x = np.arange(1,len(energies)+1)
        
        if zoom:
            energies = energies[-20:]
            x = x[-20:]

        plt.cld()
        plt.plot(x, energies, color=(215 if not neb else 37))
        plt.xlabel("Iteration #")
        plt.ylabel("Energy (kcal/mol)" if not neb else "ΔEE‡ TS(kcal/mol)")
        plt.show()
        print("\n")
        os.system(f"grep HURRAY {rootname}.out | tail -1")
        os.system(f"grep \"FINAL SINGLE POINT ENERGY\" {rootname}.out | tail -1")
        os.system(f"grep \"Total Dipole Moment\" {rootname}.out | tail -1")
        os.system(f"grep \"Magnitude\" {rootname}.out | tail -1")
        os.system(f"grep Gibbs {propname}| tail -1")
        os.system(f"grep \"G-E(el)\" {rootname}.out | tail -1")

        if scan:
            last_geom = getoutput(f"grep \"Storing optimized geometry in\" {rootname}.out | tail -1 | grep \"[0-9]\+\" -1 -o").split("\n")[-1] or 0
            
            try:
                total = getoutput(f"grep \"B [0-9]\+ [0-9]\+\" {rootname}.out").split()[-1]
            except IndexError:
                total = getoutput(f"grep \"B [0-9]\+ [0-9]\+\" {rootname}.out")

            print(f"\nLast optimized step is {int(last_geom)}/{total}")

        if neb:
            os.system(f"grep \"Starting iterations:\" {rootname}.out -A 500 ")

        else:
            homo = getoutput(f"egrep \"^ *[0-9]* +2.0000 +-[0-9].[0-9]* +[-]*[0-9]*.[0-9]*\" {rootname}.out | tail -1 | grep -o \"\-*[0-9]\+.[0-9]\+\" | tail -1")
            lumo = getoutput(f"egrep \"^ *[0-9]* +2.0000 +-[0-9].[0-9]* +[-]*[0-9]*.[0-9]*\" {rootname}.out -A 1 | tail -1 | grep -o \"\-*[0-9]\+.[0-9]\+\" | tail -1")
            
            # if homo == "" and lumo == "":
            #     homo = getoutput(f"egrep \"^ *[0-9]* +1.0000 +-[0-9].[0-9]* +[-]*[0-9]*.[0-9]*\" {rootname}.out | tail -1 | grep -o \"\-*[0-9]\+.[0-9]\+\" | tail -1")
        
            print(f"\nHOMO: {homo} eV")
            print(f"LUMO: {lumo} eV")
            print("\n"+"_"*100+"\n")

        os.system(f"grep \"TOTAL RUN TIME\" {rootname}.out | tail -1")

    out = ""
    with open(propname, 'r') as f:
        line = f.readline()
        for n in range(100000):
            line = f.readline()
            
            if "Vibrational frequencies" in line:
                out += "\nVibrational Frequencies (first 10)\n"
                while True:
                    line = f.readline()
                    if len(line.split()) > 1:
                        out += line
                        if line.split()[0] == "10":
                            break

            if not line:
                if not out:
                    out = f"\n\nNo frequencies found in output file ({propname})."
                break
            
    print(out)

scanname = f"{rootname}.relaxscanact.dat"

if scanname in files:
    
    with open(scanname, 'r') as f:
        distances, energies = [], []
        while True:
            line = f.readline()
            if not line:
                break
            d, e = line.split()
            distances.append(float(d))
            energies.append(float(e))

        energies = np.array(energies)
        energies -= np.min(energies)
        energies *= 627.5096080305927

        plt.cld()
        plt.plot(distances, energies)
        plt.scatter(distances, energies, color='red+')

        if len(scan_indices) == 2:
            plt.xlabel(f"{scan_indices} distance (A)")

        else:
            plt.xlabel(f"{scan_indices} dihedral Angle (Degrees)")

        plt.ylabel("Energy (kcal/mol)")
        plt.show()
        print("\n"+"_"*100+"\n")

    filename = rootname + ".allxyz"
    if filename in files:

        with open(filename, "r") as f:
            lines = f.readlines()

        lines = [line.replace(">", "") for line in lines]

        filename = rootname + ".all.xyz"
        with open(filename, "w") as f:
            lines = f.writelines(lines)

        mol = ccread(filename)
        maximum_id = np.argmax(energies)

        with open(f"{rootname}_scan_max.xyz", "w") as f:
            if len(scan_indices) == 2:
                units = 'dA'
            else:
                units = 'θ°'
            write_xyz(mol.atomcoords[maximum_id], mol.atomnos, f, f"Scan Maximum, {units[0]} = {round(distances[maximum_id], 3)} {units[1]}")

        print(f"Extracted scan maximum ({units[0]} = {round(distances[maximum_id], 3)} {units[1]}) to {rootname}_scan_max.xyz")
        print(f"Barrier height is {round(np.max(energies)-np.min(energies), 2)} kcal/mol")
        print("\n")

    else:
        print(f"No {filename} found in the current folder.")