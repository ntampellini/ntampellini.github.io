####################################################
options = {
    "csearch": False,
    "opt" : True,

    "method" : [#"GFN-FF",
                "AIMNet2/wB97M-V",
               ],
    "solvent" : "chloroform",
    "charge" : 0,

    # specifies an atomic index to be forced planar
    "planar" : None,
    "constrain_string" : None,
}
####################################################

import os
import sys
from time import perf_counter

import numpy as np
from firecode.algebra import dihedral, norm_of, point_angle
from firecode.calculators._xtb import xtb_gsolv, xtb_opt
from firecode.utils import Constraint, read_xyz, time_to_string, write_xyz
from InquirerPy import inquirer
from InquirerPy.base.control import Choice

from utils import multiplicity_check, get_ts_d_estimate

if len(sys.argv) == 1:
    print("\n  Optimizes the specified geometry/ies, compares results and replaces the input file. Syntax:\n\n" +
          "  python xtbopt.py filename*.xyz [newfile] [ff] [sp] [c[=file.txt]] [charge=n] [flat] [flat=]\n\n" + 
          "  filename*.xyz: Base name of input geometry file(s)\n" +
          "  newfile: Optional, creates a new file, preserving the original\n" +
          "  ff: Optional, use GFN-FF instead of GFN2-XTB\n" +
          "  sp: Optional, do not optimize but just run a single point energy calculation\n" +
          "  c: Optional, specify one or more distance/dihedral constraints to be imposed, or read from file\n" +
          "  charge=n: Optional, where \"n\" is an integer. Specify the charge to be passed to XTB\n"
        )
    quit()

allchars = ''.join(sys.argv[1:])
auto_charges = False
if '+' in allchars or '-' in allchars:
    auto_charges = inquirer.confirm(
        message="Found charge signs in filenames (+/-). Auto assign charges in input files?",
        default=True,
    ).execute()

if auto_charges:
    options["charge"] = 'auto'

options["method"] = [inquirer.select(
        message="Which level of theory would you like to use?:",
        choices=(
            Choice(value='AIMNet2/wB97M-V', name='AIMNet2/wB97M-V'),
            Choice(value='GFN2-xTB', name='GFN2-xTB (XTB)'),
            Choice(value='GFN-FF', name='GFN-FF (XTB)'),
            Choice(value='GFN2-xTB (TBLITE)', name='GFN2-xTB (TBLITE)'),
        ),
        default=options['method'],
    ).execute()]

calc = 'XTB'

if 'AIMNet2/wB97M-V' in options["method"]:
    from aimnet2_firecode.interface import aimnet2_opt, get_aimnet2_calc
    print("--> using aimnet2 via ASE")
    ase_calc = get_aimnet2_calc()
    calc = 'AIMNET2'

if 'GFN2-XTB (TBLITE)' in options["method"]:
    from firecode.ase_manipulations import ase_tblite_opt, get_ase_calc
    print("--> using tblite via ASE")
    ase_calc = get_ase_calc(('TBLITE', 'GFN2-xTB', None, options['solvent']))
    calc = 'TBLITE'

constraints, distances = [], []
newfile = False
add_free_energy = False
aimnet2 = False
smarts_string = None

if "solvent" in [kw.split('=')[0] for kw in sys.argv]:

    options["solvent"] = next(kw.split('=')[1] for kw in sys.argv if kw.split('=')[0] == 'solvent')
    print(f"--> Setting solvent to {options['solvent']}")
    sys.argv.remove(f"solvent={options['solvent']}")

if "flat" in sys.argv:
    sys.argv.remove("flat")
    assert len(sys.argv) == 2, "Single molecule only!"

    from firecode.utils import graphize
    from networkx import cycle_basis

    from utils import cycle_to_dihedrals, get_exocyclic_dihedrals

    mol = read_xyz(sys.argv[1])
    graph = graphize(mol.atomcoords[0], mol.atomnos)

    cycles = [l_ for l_ in cycle_basis(graph) if len(l_) in (7, 8, 9)]
    # if len(cycles) == 1:
    assert len(cycles) == 1, "Only 7/8/9-membered ring flips at the moment"

    print(f"--> flat: Found {len(cycles[0])}-membered ring")
    dihedrals = cycle_to_dihedrals(cycles[0])
    exocyclic = get_exocyclic_dihedrals(graph, cycles[0])
    target_angles = np.array([0 for _ in dihedrals] + [180 for _ in exocyclic])

    for (a, b, c, d), target in zip((dihedrals + exocyclic), target_angles):
        constraints.append([a, b, c, d, target])

if "flat" in [kw.split("=")[0] for kw in sys.argv]:

    flat_index = int(next(kw.split('=')[1] for kw in sys.argv if kw.split('=')[0] == 'flat'))
    print(f"--> Flattening atom {flat_index}")
    sys.argv.remove(f"flat={flat_index}")

    from firecode.graph_manipulations import neighbors
    from firecode.utils import graphize

    mol = read_xyz(sys.argv[1])
    graph = graphize(mol.atomcoords[0], mol.atomnos)
    nb = neighbors(graph, flat_index)

    assert len(nb) == 3, f'--> index {flat_index} has {len(nb)} neighbors - can only flatten trigonal pyramidal atoms.'
    
    a, b, c = nb
    constraints.append([a, flat_index, b, c, 180])
    constraints.append([b, flat_index, c, a, 180])
    constraints.append([c, flat_index, b, a, 180])
    
if "c" in sys.argv:
    
    sys.argv.remove("c")

    while True:

        data = input("Constrained indices [enter to stop]: ").split()

        if not data:
            break

        elif data[-1] == "ts":
            data[-1] = str(get_ts_d_estimate(sys.argv[1], (int(i) for i in data[0:2])))

        # try:

        assert len(data) in (2, 3, 4, 5), "Only 2-4 indices as ints + optional target as a float"

        value = None
        if '.' in data[-1]:
            value = float(data.pop(-1))

        constraint = Constraint([int(i) for i in data], value=value)                
        constraints.append(constraint)
            
        # except Exception as e:
        #     print(e)

    print(f"Specified {len(constraints)} constraints")

# constraint file
if "c" in [kw.split("=")[0] for kw in sys.argv]:
    c_filename = next((kw.split("=")[-1] for kw in sys.argv if "c=" in kw))
    sys.argv.remove(f"c={c_filename}")

    # set constraints from file
    with open(c_filename, 'r') as f:
        lines = f.readlines()

    # see if we are pattern matching
    if lines[0].startswith('SMARTS'):
        smarts_string = lines.pop(0).lstrip('SMARTS ')
        print(f'--> SMARTS line found: will pattern match and interpret constrained indices relative to the pattern')

    for line in lines:
        data = line.split()
        try:

            assert len(data) in (2, 3, 4, 5), "Only 2-4 indices as ints + optional target as a float"

            value = None
            if '.' in data[-1]:
                value = float(data.pop(-1))

            constraint = Constraint([int(i) for i in data], value=value)                
            constraints.append(constraint)

        except Exception as e:
            print(e)

    print(f'--> Read constraints from {c_filename}')

if "g" in sys.argv:
    sys.argv.remove("g")
    add_free_energy = True
    from firecode.calculators._xtb import xtb_get_free_energy
    print('--> Requested free energy calculation - adding Gcorr from GFN-FF')

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

energies, names_confs = [], []

for i, name in enumerate(names):
    try:
        data = read_xyz(name)

        print()

        outname = name if not newfile else name[:-4] + "_xtbopt.xyz"
        if newfile and (outname in os.listdir()):
            os.remove(outname)
        write_type = 'a' if newfile else 'w'

        # set charge
        if auto_charges:
            charge = name.count("+") - name.count("-")
        else:
            charge = options['charge']

        # set multiplicity
        if multiplicity_check(name[:-4], int(charge)):
            mult = 1
        else:
            mult = inquirer.text(
                message=f'It appears {name} is not a singlet. Please specify multiplicity:',
                validate=lambda inp: inp.isdigit() and int(inp) > 1,
                default="2",
            ).execute()

        for c_n, coords in enumerate(data.atomcoords):

            constrained_indices, constrained_distances, constrained_dihedrals, constrained_dih_angles = [], [], [], []
            constrained_angles_indices, constrained_angles_values = [], []

            for constraint in constraints:

                if smarts_string is not None:
                    # save original indices to revert them later, for the next conformer/molecule
                    constraint.old_indices = constraint.indices[:]

                    # correct indices from relative to the SMARTS
                    # string to absolute for this molecule
                    constraint.convert_constraint_with_smarts(coords, data.atomnos, smarts_string)
                    
                if constraint.type == 'B':

                    a, b = constraint.indices
                    if constraint.value is None:
                        constraint.value = norm_of(coords[a]-coords[b])

                    constrained_indices.append(constraint.indices)
                    constrained_distances.append(constraint.value)

                    print(f"CONSTRAIN -> d({a}-{b}) = {round(norm_of(coords[a]-coords[b]), 3)} A at start of optimization (target is {round(constraint.value, 3)} A)")

                elif constraint.type == 'A':

                    a, b, c = constraint.indices
                    if constraint.value is None:
                        constraint.value = point_angle(coords[a],coords[b],coords[c])

                    constrained_angles_indices.append(constraint.indices)
                    constrained_angles_values.append(constraint.value)
                    
                    print(f"CONSTRAIN ANGLE -> Angle({a}-{b}-{c}) = {round(point_angle(coords[a],coords[b],coords[c]), 3)}° at start of optimization, target {round(constraint.value, 3)}°")

                elif constraint.type == 'D':
                    
                    a, b, c, d = constraint.indices
                    if constraint.value is None:
                        constraint.value = dihedral(np.array([coords[a],coords[b],coords[c], coords[d]]))

                    constrained_dihedrals.append(constraint.indices)
                    constrained_dih_angles.append(constraint.value)

                    print(f"CONSTRAIN DIHEDRAL -> Dih({a}-{b}-{c}-{d}) = {round(dihedral(np.array([coords[a],coords[b],coords[c], coords[d]])), 3)}° at start of optimization, target {round(constraint.value, 3)}°")
            
            action = "Optimizing" if options["opt"] else "Calculating SP energy on"
            for method in options["method"]:
                print(f'{action} {name} - {i+1} of {len(names)}, conf {c_n+1} of {len(data.atomcoords)} ({method}) - CHG={charge} MULT={mult}')
                t_start = perf_counter()

                if calc == 'AIMNET2':

                    if mult != 1:
                        raise NotImplementedError('AIMNET2 only supports singlets.')

                    coords, energy, success = aimnet2_opt(
                                                            coords,
                                                            data.atomnos,

                                                            constrained_indices=constrained_indices,
                                                            constrained_distances=constrained_distances,

                                                            constrained_angles_indices=constrained_angles_indices,
                                                            constrained_angles_values=constrained_angles_values,

                                                            constrained_dihedrals_indices=constrained_dihedrals,
                                                            constrained_dihedrals_values=constrained_dih_angles,

                                                            ase_calc=ase_calc,
                                                            traj=name[:-4] + "_traj",
                                                            logfunction=print,
                                                            charge=charge,
                                                            maxiter=500 if options["opt"] else 1,
                                                            solvent=options["solvent"],

                                                            # title='XTBOPT_temp',
                                                            # debug=True,
                                                        )
                    
                elif calc == "XTB":
                
                    coords, energy, success = xtb_opt(
                                                    coords,
                                                    data.atomnos,
                                                    constrained_indices=constrained_indices,
                                                    constrained_distances=constrained_distances,

                                                    constrained_dihedrals=constrained_dihedrals,
                                                    constrained_dih_angles=constrained_dih_angles,

                                                    constrained_angles_indices=constrained_angles_indices,
                                                    constrained_angles_values=constrained_angles_values,

                                                    constrain_string=options["constrain_string"],
                                                    method=method,
                                                    solvent=options["solvent"],
                                                    charge=charge,
                                                    mult=mult,
                                                    opt=options["opt"],
                                                    )

                elif calc == 'TBLITE':

                    coords, energy, success = ase_tblite_opt(
                                                    coords,
                                                    data.atomnos,
                                                    tb_calc=ase_calc,
                                                    method=method,

                                                    constrained_indices=constrained_indices,
                                                    constrained_distances=constrained_distances,

                                                    constrained_angles_indices=constrained_angles_indices,
                                                    constrained_angles_values=constrained_angles_values,

                                                    constrained_dihedrals_indices=constrained_dihedrals,
                                                    constrained_dihedrals_values=constrained_dih_angles,

                                                    charge=charge,
                                                    mult=1,
                                                    solvent=options["solvent"],
                                                    
                                                    maxiter=500 if options["opt"] else 1,
                                                    traj=name[:-4] + "_traj",
                                                    logfunction=print,

                                                    # title='temp',
                                                    # debug=False,
                                                )

                elapsed = perf_counter() - t_start

                if energy is None:
                    print(f'--> ERROR: Optimization of {name} crashed. ({time_to_string(elapsed)})')

                elif options["opt"]:
                    with open(outname, write_type) as f:
                        write_xyz(coords, data.atomnos, f, title=f'Energy = {energy} kcal/mol')
                    print(f"{'Appended' if write_type == 'a' else 'Wrote'} optimized structure at {outname} - {time_to_string(elapsed)}\n")

            if add_free_energy:
                sph = (len(constraints) != 0)
                print(f'Calculating Free Energy contribution{" (SPH)" if sph else ""} on {name} - {i+1} of {len(names)}, conf {c_n+1} of {len(data.atomcoords)} ({method})')
                gcorr = xtb_get_free_energy(coords, data.atomnos, method='GFN-FF', solvent=options["solvent"], charge=options["charge"], sph=sph, grep='Gcorr')
                print(f'GCORR: {name}, conf {c_n+1} - {gcorr:.2f} kcal/mol')
                energy += gcorr

            energies.append(energy)
            names_confs.append(name[:-4]+f"_conf{c_n+1}")

    except RuntimeError as e:
        print("--> ", name, " - ", e)
        continue

    if constraints:
        print(f'Constraints: final values')

        for constraint in constraints:
            if constraint.type == 'B':
                a, b = constraint.indices
                final_value = norm_of(coords[a]-coords[b])
                uom = ' Å'
            
            elif constraint.type == 'A':
                a, b, c = constraint.indices
                final_value = point_angle(coords[a],coords[b],coords[c])
                uom = '°'

            elif constraint.type == 'D':
                a, b, c, d = constraint.indices
                final_value = dihedral(np.array([coords[a],coords[b],coords[c], coords[d]]))
                uom = '°'

            indices_string = '-'.join([str(i) for i in constraint.indices])
            print(f"CONSTRAIN -> {constraint.type}({indices_string}) = {round(final_value, 3)}{uom}")

            # revert original indices for the next molecule
            if smarts_string is not None:
                constraint.indices = constraint.old_indices

        print()
if None not in energies:

    if len(names_confs) > 1:
        min_e = min(energies)
    else:
        min_e = 0

    print()
    longest_name = max([len(s) for s in names_confs])
    for nc, energy in zip(names_confs, energies):
        print(f"ENERGY -> {nc:{longest_name}s} = {energy-min_e:.2f} kcal/mol")