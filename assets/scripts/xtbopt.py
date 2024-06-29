####################################################
options = {
    "csearch": False,
    "opt" : True,

    "method" : [#"GFN-FF",
                "GFN2-XTB",
               ],
    "solvent" : "ch2cl2",
    "charge" : 0,

    # specifies an atomic index to be forced planar
    "planar" : None,
    "constrain_string" : None,
}
####################################################

from firecode.calculators._xtb import xtb_opt, xtb_gsolv
from firecode.utils import read_xyz, write_xyz, time_to_string
from firecode.algebra import norm_of, dihedral
import os 
import numpy as np
import sys
from time import perf_counter

if len(sys.argv) == 1:
    print("\n  Optimizes the specified geometry/ies, compares results and replaces the input file. Syntax:\n\n" +
          "  python xtbopt.py filename*.xyz [newfile] [ff] [sp] [c] [charge=n]\n\n" + 
          "  filename*.xyz: Base name of input geometry file(s)\n" +
          "  newfile: Optional, creates a new file, preserving the original\n" +
          "  ff: Optional, use GFN-FF instead of GFN2-XTB\n" +
          "  sp: Optional, do not optimize but just run a single point energy calculation\n" +
          "  c: Optional, specify one or more distance/dihedral constraints to be imposed\n" +
          "  charge=n: Optional, where \"n\" is an integer. Specify the charge to be passed to XTB\n"
        )
    quit()

if "ff" in sys.argv:
    sys.argv.remove("ff")
    options["method"] = ["GFN-FF"]
    print("--> Using GFN-FF force field")

constraints, distances = [], []
newfile = False
free_energy = False
aimnet2 = False

if "solvent" in [kw.split('=')[0] for kw in sys.argv]:

    options["solvent"] = next(kw.split('=')[1] for kw in sys.argv if kw.split('=')[0] == 'solvent')
    print(f"--> Setting solvent to {options['solvent']}")
    sys.argv.remove(f"solvent={options['solvent']}")

if "aim" in sys.argv:

    from aimnet2_firecode.interface import get_aimnet2_calc, aimnet2_opt
    
    print("--> using aimnet2 via ASE")
    calc = get_aimnet2_calc()
    aimnet2 = True
    options["method"] = ['AIMNet2/wB97M-V']
    # options["solvent"] = None
    sys.argv.remove("aim")

if "flat" in sys.argv:
    sys.argv.remove("flat")
    assert len(sys.argv) == 2, "Single molecule only!"

    from firecode.utils import graphize
    from networkx import cycle_basis
    from utils import cycle_to_dihedrals, get_exocyclic_dihedrals

    mol = read_xyz(sys.argv[1])
    graph = graphize(mol.atomcoords[0], mol.atomnos)

    cycles = [l_ for l_ in cycle_basis(graph) if len(l_) in (7, 8, 9)]
    assert len(cycles) == 1, "Only 7/8/9-membered ring flips at the moment"

    dihedrals = cycle_to_dihedrals(cycles[0])
    exocyclic = get_exocyclic_dihedrals(graph, cycles[0])
    target_angles = np.array([0 for _ in dihedrals] + [180 for _ in exocyclic])

    for (a, b, c, d), target in zip((dihedrals + exocyclic), target_angles):
        constraints.append([a, b, c, d, target])

    print(f"--> flat: Found {len(cycles[0])}-membered ring")

if "c" in sys.argv:
    sys.argv.remove("c")

    while True:

        data = input("Constrained indices [enter to stop]: ").split()

        if data == []:
            break

        try:

            assert len(data) in (2, 3, 4, 5), "Only distances (+optional target, 2 or 3 numbers) or dihedrals (+optional target) supported (4 or 5 numbers)"

            if len(data) == 2:
                indices = [int(x) for x in data]
            
            elif len(data) == 3:
                indices = [int(x) for x in data[0:2]]+ [float(data[2])]

            elif len(data) == 4:
                indices = [int(x) for x in data]

            elif len(data) == 5:
                indices = [int(x) for x in data[0:4]] + [float(data[4])]
                
            constraints.append(indices)
            
        except Exception as e:
            print(e)

    print(f"Specified {len(constraints)} constraints")

if "g" in sys.argv:
    sys.argv.remove("g")
    free_energy = True
    from firecode.calculators._xtb import xtb_get_free_energy
    print('--> Requested free energy calculation')


if any(["constr" in kw.split(",") for kw in sys.argv]):
    string = [kw for kw in sys.argv if "constr" in kw.split(",")][0]
    sys.argv.remove(string)
    # clunky way to do it, I know ... but allows for non-interactive usage. "constr,i1,i2,ds"
    
    parts = string.split(",")
    assert len(parts) == 4
    distances.append(float(parts[3]))
    constraints.append([int(x) for x in parts[1:3]]+[float(parts[3])])

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

# if "planar" in sys.argv:
#     sys.argv.remove(f"planar")
#     options["constrain_string"] = "dihedral: 7, 8, 9, 15, 180\ndihedral: 9, 8, 15, 7, 180\ndihedral: 15, 8, 7, 9, 180\n force constant=1.0"
#     print("--> !PLANAR")

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

energies, names_confs, constrained_indices, constrained_distances, constrained_dihedrals, constrained_dih_angles = [], [], [], [], [], []

for i, name in enumerate(names):
    data = read_xyz(name)

    print()

    outname = name if not newfile else name[:-4] + "_xtbopt.xyz"
    if newfile and (outname in os.listdir()):
        os.remove(outname)
    write_type = 'a' if newfile else 'w'

    for c_n, coords in enumerate(data.atomcoords):
        for constraint in constraints:
            if len(constraint) in (2, 3):
                if len(constraint) == 2:
                    a, b = constraint
                    target = norm_of(coords[a]-coords[b])
                else:
                    a, b, target = constraint

                constrained_indices.append([a, b])
                constrained_distances.append(target)

                print(f"CONSTRAIN -> d({a}-{b}) = {round(norm_of(coords[a]-coords[b]), 3)} A at start of optimization (target is {target} A)")
                       
            elif len(constraint) in (4, 5):
                if len(constraint) == 4:
                    a, b, c, d = constraint
                    dih_angle = dihedral(np.array([coords[a],coords[b],coords[c], coords[d]]))
                else:
                    a, b, c, d, dih_angle = constraint

                constrained_dihedrals.append([a, b, c, d])
                constrained_dih_angles.append(dih_angle)

                print(f"CONSTRAIN DIHEDRAL -> Dih({a}-{b}-{c}-{d}) = {round(dihedral(np.array([coords[a],coords[b],coords[c], coords[d]])), 3)} ° at start of optimization, target {round(dih_angle, 3)}")
        
        action = "Optimizing" if options["opt"] else "Calculating SP energy on"
        for method in options["method"]:
            print(f'{action} {name} - {i+1} of {len(names)}, conf {c_n+1} of {len(data.atomcoords)} ({method})')
            t_start = perf_counter()

            if aimnet2:

                assert constrained_dihedrals == [], "No dih angles for aim yet..."

                coords, energy, success = aimnet2_opt(
                                                        coords,
                                                        data.atomnos,
                                                        constrained_indices=constrained_indices,
                                                        constrained_distances=constrained_distances,
                                                        ase_calc=calc,
                                                        # traj=name[:-4] + "_traj.xyz",
                                                        logfunction=print,
                                                        title='temp',
                                                        charge=options["charge"],
                                                        maxiter=500 if options["opt"] else 1,
                                                    )
                
                if options["solvent"] is not None:
                    gsolv = xtb_gsolv(
                                        coords,
                                        data.atomnos,
                                        charge=options['charge'],
                                        solvent=options['solvent'],
                                    )
                    print(f'--> {name}: ALPB GSolv = {gsolv:.2f} kcal/mol')
                    energy += gsolv

            else:
            
                coords, energy, success = xtb_opt(coords,
                                                data.atomnos,
                                                constrained_indices=constrained_indices,
                                                constrained_distances=constrained_distances,
                                                constrained_dihedrals=constrained_dihedrals,
                                                constrained_dih_angles=constrained_dih_angles,
                                                constrain_string=options["constrain_string"],
                                                method=method,
                                                solvent=options["solvent"],
                                                charge=options["charge"],
                                                opt=options["opt"],
                                                )

            elapsed = perf_counter() - t_start

            if energy is None:
                print(f'--> ERROR: Optimization of {name} crashed. ({time_to_string(elapsed)})')

            elif options["opt"]:
                with open(outname, write_type) as f:
                    write_xyz(coords, data.atomnos, f, title=f'Energy = {energy} kcal/mol')
                print(f"{'Appended' if write_type == 'a' else 'Wrote'} optimized structure at {outname} - {time_to_string(elapsed)}")

        if free_energy:
            sph = (len(constraints) != 0)
            print(f'Calculating Free Energy{" (SPH)" if sph else ""} on {name} - {i+1} of {len(names)}, conf {c_n+1} of {len(data.atomcoords)} ({method})')
            energy = xtb_get_free_energy(coords, data.atomnos, method=method, solvent=options["solvent"], charge=options["charge"], sph=sph)

        energies.append(energy)
        names_confs.append(name[:-4]+f"_conf{c_n+1}")

    for constraint, target in zip(constraints, distances):
        if len(constraint) == 2:
            print(f"CONSTRAIN -> d({a}-{b}) = {round(norm_of(coords[a]-coords[b]), 3)} A")

if None not in energies:

    if len(names_confs) > 1:
        min_e = min(energies)
    else:
        min_e = 0

    print()
    longest_name = max([len(s) for s in names_confs])
    for nc, energy in zip(names_confs, energies):
        print(f"ENERGY -> {nc:{longest_name}s} = {energy-min_e:.2f} kcal/mol")