import os
from subprocess import getoutput
from cclib.io import ccread
import pickle
from rich.traceback import install
install(show_locals=True)

outname = input("--> ORCA NMR extractor - pick a name : ").rstrip() + ".pickle"

startdir = os.getcwd()
all_files_present = True
all_data_present = True

# get names and create storage dictionary
data = dict()
for name in os.listdir():
    if not name.startswith("slurm") and name.endswith(".out"):
        data[name] = dict()

######################################################################################### CHECKING SP

# check that all the files are present before collecting data
n_files = len(data.keys())
print(f"Found {n_files} opt/freq files. Checking single points...")

# checking single point calcs
os.chdir("SP")
missing_files = list(data.keys())
for name in os.listdir():
    if not name.startswith("slurm") and name.endswith(".out"):
        if name in missing_files:
            missing_files.remove(name)

# printing out missing sp files
if missing_files != []:
    print(f"Error: only found {int(n_files-len(missing_files))}/{n_files} SP files!")
    for name in missing_files:
        print(f"SP: missing {name}")
    all_files_present = False
else:
    print("All single point output files found.")

# checking single point calcs
os.chdir("NMR")
missing_files = list(data.keys())
for name in os.listdir():
    if not name.startswith("slurm") and name.endswith(".out"):
        if name in missing_files:
            missing_files.remove(name)

# printing out missing sp files
if missing_files != []:
    print(f"Error: only found {int(n_files-len(missing_files))}/{n_files} NMR files!")
    for name in missing_files:
        print(f"NMR: missing {name}")
    all_files_present = False
else:
    print("All NMR output files found.")

######################################################################################### COLLECTING DATA

if all_files_present:

    print("Collecting data...")

    #get back at the start directory
    os.chdir(startdir)

    # start collecting: extract free energy correction and add it to the dict
    for name in data.keys():
        corr = float(getoutput(f"grep \"G-E(el)\" {name} | tail -1").split()[2])
        # print(f"corr - {name} - {corr}")
        data[name]["gcorr"] = corr

    # move to SP folder and extract sp energy and geometry
    os.chdir("SP")
    # sp_energies, geometries = [], []
    for name in data.keys():
        energy = float(getoutput(f"grep \"FINAL SINGLE POINT\" {name} | tail -1").split()[4])
        data[name]["sp_energy"] = energy
        data[name]["free_energy"] = energy + data[name]["gcorr"]
        # print(f"energy - {name} - {energy}")

        mol = ccread(name[:-4]+".xyz")
        data[name]["mol"] = mol
        # print(f"{name} - {len(mol.atomcoords[0])} atoms")

    os.chdir("NMR")
    for name in data.keys():
        lines = getoutput(f"grep \"CHEMICAL SHIELDING SUMMARY\" {name} -A {len(data[name]['mol'].atomcoords[0])+10}").splitlines()
        shieldings = [float(line.split()[2]) for line in lines[6:-4] if line != '']
        data[name]["shieldings"] = shieldings

######################################################################################### INTEGRITY CHECK

    n_corr = [data[name].get('gcorr', None) is not None for name in data.keys()].count(True)
    if n_corr < n_files:
        print(f"--> Only found {n_corr}/{n_files} G-E(el) corrections! Have all the calculations completed successfully?")
        all_data_present = False

        for name in data.keys():
            if type(data[name]["gcorr"]) is not float:
                print(f"OPT/FREQ data: missing G-E(el) correction from {name}")
    else:
        print("Found all G-E(el) corrections.")


    n_sp = [data[name].get('sp_energy', None) is not None for name in data.keys()].count(True)
    if n_sp < n_files:
        print(f"--> Only found {n_sp}/{n_files} SP energies! Have all the calculations completed successfully?")
        all_data_present = False

        for name in data.keys():
            if type(data[name]["sp_energy"]) is not float:
                print(f"SP data: missing single point energy from {name}")
    else:
        print("Found all single point energies.")

    n_shieldings = [data[name]['shieldings'] != [] for name in data.keys()].count(True)
    if n_shieldings < n_files:
        print(f"--> Only found {n_shieldings}/{n_files} NMR shieldings! Have all the calculations completed successfully?")
        all_data_present = False

        for name in data.keys():
            if data[name]["shieldings"] == []:
                print(f"NMR data: missing shieldings from {name}")
    else:
        print("Found all NMR shieldings.")

######################################################################################### DATA DUMP

    if all_data_present:

        # back to the place we ran the script
        os.chdir(startdir)

        with open(outname, "wb") as _f:
            pickle.dump(data, _f)
            
        print(f'--> Written {outname} ({round(os.path.getsize(outname)/1000, 1)} kB)')

    else:
        print(f"Missing some data: cannot export a complete .pickle file.")

else:
    print(f"Missing .out files: cannot export data.")