import os
import sys
from subprocess import getoutput

import numpy as np
import plotext as plt
from cclib.io import ccread
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from numpy.linalg import norm

from utils import dihedral, write_xyz

EH_TO_KCAL = 627.5096080305927

def main(loop=False):

    plt.theme("pro")
    plt.plotsize(100,25)
    plot_distance = False
    scan = False
    neb = False
    compound = False
    indices = None
    zoom = False
    scratchdir = '/vast/palmer/scratch/miller/nt383'

    # if no outname is provided, list running jobs
    if len(sys.argv) == 1 or loop:

        cwd = os.getcwd()
        scratchnames = getoutput(f'ls {scratchdir}').split()
        choices = []
        
        pids = [line.split()[0] for line in getoutput("squeue --me").split('\n')[1:]]
        for pid in pids:
            try:
                workdir = [line.split('=')[1] for line in getoutput(f"scontrol show job {pid}").split("\n") if 'WorkDir' in line][0]
                # if workdir == cwd:
                jobdir = next(name for name in scratchnames if pid in name)
                jobname = ('-'.join(jobdir.split('-')[1:]))[:-5].strip(pid).strip('_')
                fullname = scratchdir+"/"+jobdir+'/'+jobname+'.out'
                choices.append(Choice(value=fullname, name=f'{pid} - {workdir}/{jobname}.out'))

            except StopIteration:
                continue

        if choices:
            loop = True
            fullname = inquirer.select(
                message="Which running job would you like to check?",
                choices=choices+[Choice(value=None, name="(quit)")],
                default=choices[0].value,
                ).execute()
            
            if fullname == None:
                sys.exit()

            sys.argv.append(fullname)
            print(f'--> Checking PID {pid} in scratch directory {fullname}')

        else:
            print("\n  Check running and completed ORCA jobs and prints/plots relevant information. " +
                "For scans, extracts local maxima structures. Syntax:\n\n" +
                "  python check.py filename[.out] [PID] [i1] [i2] [zoom]\n\n" + 
                "  filename/PID: base name of input/output file (with or without extension),\n" +
                "    or alternatively, the job PID (needs the scratch directory to be set)\n" +
                "  i1/i2: optional, plot i1-i2 distance during optimization steps\n" +
                "  zoom: optional, only shows plot of last 20 iterations.\n")
            quit()

    def is_nonempty(string) -> bool:
        '''
        returns True if a string contains more
        than spaces and newline chars.

        '''
        return not not string.replace('\n', '').strip()

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

    if rootname.endswith(".xyz"):

        energies = [float(line.split()[-1]) for line in getoutput(f'grep \'E.*-\d*.\' {rootname}').split('\n')]
        energies_kcal = np.array(energies) * EH_TO_KCAL
        energies_kcal -= min(energies_kcal)

        plt.cld()
        plt.plot(energies_kcal, color=37)
        plt.scatter(energies_kcal, color='red+')
        plt.xlabel("Structure #")
        plt.ylabel("Rel. Energy (kcal/mol)")
        plt.show()

        print(f'Min: Structure {energies.index(min(energies))}/{len(energies)} ({min(energies):.8f} Eh)')
        print(f'Max: Structure {energies.index(max(energies))}/{len(energies)} ({min(energies):.8f} Eh)')

        return
    else:
        rootname = rootname.split(".")[0]

    # catch use of PID to check a job
    if rootname.isdigit():
        scratch_folder = getoutput(f'ls {scratchdir} | grep {rootname}').rstrip('\n')

        if scratch_folder == '':
            print(f'Job {rootname} not found in scratch.')
            sys.exit()

        os.chdir(os.path.join(scratchdir, scratch_folder))
        print(f'--> Checking PID {rootname} from scratch directory {os.getcwd()}')
        rootname = getoutput(f'ls *.out').rstrip('\n').split('.')[0]

    inpname = f"{rootname}.inp"
    files = os.listdir()

    if inpname in files:
        with open(inpname, "r") as f:
            while True:

                line = f.readline()
                if not line:
                    break

                if "%GEOM SCAN" in line.upper():
                    scan = True
                    frags = f.readline().upper().split()

                    if frags[0] == "B":
                        scantype = "distance"
                        scan_indices = tuple([int(i) for i in frags[1:3]])

                    elif frags[0] == "A":
                        scantype = "angle"
                        scan_indices = tuple([int(i) for i in frags[1:4]])

                    elif frags[0] == "D":
                        scantype = "dihedral"
                        scan_indices = tuple([int(i) for i in frags[1:5]])

                    else:
                        scan = False

                if "NEB" in line.upper():
                    neb = True

                if ".CMP" in line.upper():
                    compound = True

    filename = rootname + "_trj.xyz"

    if filename in files and indices is not None:

        mol = ccread(filename)

        if len(indices) == 2:
            y = [norm(c[indices[0]]-c[indices[1]]) for c in mol.atomcoords]
            tag = "distance (Å)"
            uom = "Å"

        elif len(indices) == 4:
            y = [dihedral(c[indices]) for c in mol.atomcoords]
            tag = "dihedral angle (Degrees)"
            uom = "°"

        else:
            raise Exception("Provide 2 or 4 indices")

        plot = plt.plot(y)
        plt.xlabel("Iteration #")
        plt.ylabel(f"{indices} {tag}")
        plt.show()
        print(f"\n Last geometry: {y[-1]:.3f} {uom}")
        print("\n"+"_"*100+"\n")

    propnames = (
        f"{rootname}.property.txt",
        f"{rootname}_Compound_1.property.txt",
        f"{rootname}_Compound_2.property.txt",
        f"{rootname}_Compound_3.property.txt",
    )

    for propname in propnames:
        if propname in files:
            
            if not neb:
                # with open(propname, 'r') as f:
                #     energies = []
                #     while True:
                #         line = f.readline()
                #         if "Total DFT Energy" in line:
                #             energies.append(float(line.split()[6]))

                #         if not line:
                #             break
                try:
                    energies = [float(line.split()[3]) for line in getoutput(f'grep \'&FINALEN\' {propname}').split('\n')]
                except IndexError:
                    energies = []

            else:
                # energies = []
                # lines = getoutput(f"grep \"Starting iterations:\" {rootname}.out -A 500 ").splitlines()
                # # energies = [float(line.split()[3]) for line in lines[4:]]
                # for line in lines[4:]:
                #     try:
                #         assert line.split()[0] != "Convergence"
                #         energies.append(float(line.split()[3]))

                #     except (IndexError, ValueError, AssertionError):
                #         continue

                raise NotImplementedError()
        
            if len(energies) > 0:
                last_E = energies[-1]
                energies = np.array(energies)
                energies -= np.min(energies) if not neb else 0
                energies *= EH_TO_KCAL
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
                os.system(f"grep \"G-E(el)\" {rootname}.out | tail -1")
                
                if (g := getoutput(f"grep \'&FREEENERGYG\' {propname}| tail -1")) != "":
                    print(f'FINAL FREE ENERGY    {float(g.split()[3])} Eh')

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

            # PRINT VIBRATIONAL FREQUENCIES

            vib_block = getoutput(f'grep \'&FREQ \' {propname} -A 20')

            if vib_block != "":

                vib_blocks = [block for block in vib_block.split('&FREQ') if is_nonempty(block)]
                for b, block in enumerate(vib_blocks):
                    try:
                        freqs = [float(line.split()[1]) for line in block.split('\n')[3:] if len(line.split()) > 1]
                        # freqs = [n for n in freqs if abs(n) > 1E-3]
                        print(f"\nVibrational Frequencies (first 10) - [{b+1}/{len(vib_blocks)}]\n")
                        for i, freq in enumerate(freqs):
                            print(f'  {i:2}    {freq:7.2f} cm^-1')
                            if i == 10:
                                break
                        print('\n')

                        # if one is a TS mode, characterize it
                        mode_block = getoutput('grep -E "^\s*[0-9]+\.\s+B\([A-Z]\s*[0-9]+,[A-Z]\s*[0-9]+\)\s+([-+]?[0-9]*\.[0-9]+)\s+([-+]?' +
                            f'[0-9]*\.[0-9]+)\s+([-+]?[0-9]*\.[0-9]+)\s+([-+]?[0-9]*\.[0-9]+)\s+([-+]?[0-9]*\.[0-9]+)\s*$" {rootname}.out')
                        
                        if mode_block != '':
                            mode_lines = [line for line in mode_block.split("\n")]
                            mode_descs = {' '.join(line.split()[1:4]) for line in mode_lines}
                            print('--> TS mode characterization')
                            for mode_desc in mode_descs:
                                # only pick the last value for that mode 
                                value = [line.split()[7] for line in mode_lines if mode_desc in line][-1]
                                print(f"  {mode_desc}    {value}")
                            print("\n")


                    except IndexError:
                        pass
            else:
                print(f'No frequencies found in output file ({propname}).\n')

            # out = ""
            # with open(propname, 'r') as f:
            #     line = f.readline()
            #     for n in range(100000):
            #         line = f.readline()
                    
            #         if "Vibrational frequencies" in line:
            #             out += "\nVibrational Frequencies (first 10)\n"
            #             while True:
            #                 line = f.readline()
            #                 if len(line.split()) > 1:
            #                     out += line
            #                     if line.split()[0] == "10":
            #                         break

            #         if not line:
            #             if not out:
            #                 out = f"\n\nNo frequencies found in output file ({propname})."
            #             break
                    
            # print(out)

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
            energies *= EH_TO_KCAL

            plt.cld()
            plt.plot(distances, energies)
            plt.scatter(distances, energies, color='red+')

            if len(scan_indices) == 2:
                plt.xlabel(f"{scan_indices} distance (Å)")

            if len(scan_indices) == 3:
                plt.xlabel(f"{scan_indices} planar angle (degrees)")

            else:
                plt.xlabel(f"{scan_indices} dihedral angle (degrees)")

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

    # else:
    #     lines = getoutput(f'grep \"Total DFT Energy\" {rootname}.property.txt')
    #     energies = [float(line.split()[6]) for line in lines]

    if compound:
        os.system(f"grep \'\-\->\' {rootname}.out")

    # if no outname was provided, loop
    if loop:
        sys.argv = [sys.argv[0]]
        main(loop=loop)

if __name__ == '__main__':
    main()