####################################################
options = {
    "csearch": False,
    "opt" : True,

    "method" : [#"GFN-FF",
                "GFN2-XTB",
               ],
    "solvent" : "ch2cl2",
    "charge" : 0,
}
####################################################

from tscode.calculators._xtb import xtb_opt
from tscode.utils import read_xyz, write_xyz
from tscode.algebra import norm_of, dihedral
import os 
import numpy as np
import sys

if len(sys.argv) == 1:
    print(f"\n  Optimizes the specified geometry/ies, compares results and replaces the input file. Syntax:\n\n" +
           "  python xtbopt.py filename*.xyz [newfile] [ff] [sp] [c] [charge]\n\n" + 
           "  filename*.xyz: base name of input geometry file(s)\n" +
           "  newfile: optional, creates a new file, preserving the original\n" +
           "  ff: optional, use GFN-FF instead of GFN2-XTB\n" +
           "  sp: optional, do not optimize but just run a single point energy calculation\n" +
           "  c: optional, specify one or more distance/dihedral constraints to be imposed\n" +
           "  charge: optional, specify the charge to be passed to XTB\n"
           )
    quit()

if "ff" in sys.argv:
    sys.argv.remove("ff")
    options["method"] = ["GFN-FF"]
    print("--> Using GFN-FF force field")

constraints, distances = [], []
newfile = False
free_energy = False

if "c" in sys.argv:
    sys.argv.remove("c")

    while True:

        data = input("Constrained indices [enter to stop]: ").split()

        if data == []:
            break

        try:

            assert len(data) in (2, 3, 4), "Only distances (+optional target) supported (2 or 3 numbers) or dihedrals (4 numbers)"

            if len(data) == 3:
                distances.append(float(data[2]))
                data = data[:2]

            else:
                distances.append(None)

            indices = [int(x) for x in data]
            constraints.append(indices)
            
        except Exception as e:
            print(e)

    print(f"Specified {len(constraints)} constraints")

if "g" in sys.argv:
    sys.argv.remove("g")
    free_energy = True
    from tscode.calculators._xtb import xtb_get_free_energy
    print('--> Requested free energy calculation')


if any(["constr" in kw.split(",") for kw in sys.argv]):
    string = [kw for kw in sys.argv if "constr" in kw.split(",")][0]
    sys.argv.remove(string)
    # clunky way to do it, I know ... but allows for non-interactive usage. "constr,i1,i2,ds"
    
    parts = string.split(",")
    assert len(parts) == 4
    distances.append(float(parts[3]))
    constraints.append([int(x) for x in parts[1:3]])

if "sp" in sys.argv:
    sys.argv.remove("sp")
    options["opt"] = False
    print("--> Single point calculation requested (no optimization)")

if "newfile" in sys.argv:
    sys.argv.remove("newfile")
    newfile = True
    print("printing to new files")

if "charge" in [kw.split("=")[0] for kw in sys.argv]:
    options["charge"] = next((kw.split("=")[-1] for kw in sys.argv if "charge" in kw))
    sys.argv.remove(f"charge={options['charge']}")

for option, value in options.items():
    print(f"--> {option} = {value}")

os.chdir(os.getcwd())
names = []

if len(sys.argv) < 2:
    for f in os.listdir():
        if os.path.isfile(f) and f.split('.')[1] == 'xyz':
            names.append(f)
else:
    for name in sys.argv[1:]:
        names.append(name)

energies, names_confs = [], []

for i, name in enumerate(names):
    data = read_xyz(name)

    outname = name if not newfile else name[:-4] + "_xtbopt.xyz"
    if newfile and (outname in os.listdir()):
        os.remove(outname)
    write_type = 'a' if newfile else 'w'

    for c, coords in enumerate(data.atomcoords):
        for constraint, target in zip(constraints, distances):
            if len(constraint) == 2:
                a, b = constraint
                print(f"CONSTRAIN -> d({a}-{b}) = {round(norm_of(coords[a]-coords[b]), 3)} A at start of optimization (target is {target} A)")
            else:
                a, b, c, d = constraint
                print(f"CONSTRAIN DIHEDRAL -> Dih({a}-{b}-{c}-{d}) = {round(dihedral(np.array([coords[a],coords[b],coords[c], coords[d]])), 3)} ° at start of optimization")
        
        action = "Optimizing" if options["opt"] else "Calculating SP energy on"
        for method in options["method"]:
            print(f'{action} {name} - {i+1} of {len(names)}, conf {c+1} of {len(data.atomcoords)} ({method})')
            coords, energy, success = xtb_opt(coords, data.atomnos, constrained_indices=constraints, constrained_distances=distances, method=method, solvent=options["solvent"], charge=options["charge"], opt=options["opt"])

            if options["opt"]:
                with open(outname, write_type) as f:
                    write_xyz(coords, data.atomnos, f, title=f'Energy = {energy} kcal/mol')
                print(f"{'Appended' if write_type == 'a' else 'Wrote'} optimized structure at {outname}")

        if free_energy:
            print(f'Calculating Free Energy on {name} - {i+1} of {len(names)}, conf {c+1} of {len(data.atomcoords)} ({method})')
            energy = xtb_get_free_energy(coords, data.atomnos, method=method, solvent=options["solvent"], charge=options["charge"])

        energies.append(energy)
        names_confs.append(name[:-4]+f"_conf{c+1}")

    for constraint, target in zip(constraints, distances):
        if len(constraint) == 2:
            print(f"CONSTRAIN -> d({a}-{b}) = {round(norm_of(coords[a]-coords[b]), 3)} A")

if len(names_confs) > 1:
    min_e = min(energies)
else:
    min_e = 0

for nc, energy in zip(names_confs, energies):
    print(f"ENERGY -> {nc} = {round(energy-min_e, 3)} kcal/mol")