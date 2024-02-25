from subprocess import getoutput
import sys

if len(sys.argv) == 1:
    print(f"\n  Compare completed ORCA jobs electronic/free energy of groups of jobs and prints " +
           "collected energy values. Syntax:\n\n" +
           "  python compare.py conf*.out [g] [x=folder/file.xyz]\n\n" + 
           "  conf*.out: base name with generic asterisk descriptor (group of output files)\n" +
           "  g: optional, uses free energy for the comparison (needs freq calculations).\n" +
           "  x: extract conformers from the specified files, cutting at 5 kcal/mol and removing duplicates.\n" +
           "     Writes structures to folder/file.xyz")
    quit()

g = False
if "g" in sys.argv:
    g = True
    look_for = "Final Gibbs free energy"
    sys.argv.remove("g")
else:
    look_for = "FINAL SINGLE POINT ENERGY"

x = False
x_thr = 5
if "x" in [kw.split("=")[0] for kw in sys.argv]:
    outname = next((kw.split("=")[-1] for kw in sys.argv if "x" in kw))
    sys.argv.remove(f"x={outname}")
    x = True

names, energies, corrections = [], [], []
for name in sys.argv[1:]:
    print(f'Reading {name}...', end="\r")
    try:
        energy = float(getoutput(f'grep \"{look_for}\" {name} | tail -1').rstrip(" Eh").split()[-1])

        if g:
            gcorr = float(getoutput(f'grep \"G-E(el)\" {name} | tail -1').split()[2])
            corrections.append(gcorr)

        names.append(name)
        energies.append(energy)
    except IndexError:
        pass

print()
min_e = min(energies)
names, energies = zip(*sorted(zip(names, energies), key=lambda x: x[1]))

print(f"\nGrepped {look_for}\n")
maxlen = max([len(name) for name in names]) + 2

for name, energy in zip(names, energies):
    converged = True if getoutput(f'grep HURRAY {name}') != '' else False
    iterations = "conv" if converged else getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {name} -c')
    print(f'{name:<{maxlen}} - {round((energy-min_e)*627.509608030593, 3): 10} - [{iterations}{"" if converged else " iterations"}]')

lowest = energies.index(min(energies))
print(f'\nBest is {names[lowest]} at {energies[lowest]}\n')

table_title = f'Name                 {"Free" if g else "Electronic"} Energies (Eh)'
if g:
    table_title += "  G-E(el) (Eh)"
print(table_title)

for i, (name, energy) in enumerate(zip(names, energies)):
    s = '{:20}    {}'.format(name, energy)

    if g:
        s += " "*3 + str(corrections[i])

    print(s)

if x:
    from tscode.utils import read_xyz, write_xyz, graphize
    from prune import cl_similarity_refining
    from update import update
    import numpy as np

    low_e_names = [name for n, name in enumerate(names) if (energies[n]-min_e)*627.509608030593 < x_thr]

    print(f"Keeping {len(low_e_names)}/{len(names)} structures ({x_thr} kcal/mol thr)")

    print(f'\nUpdating optimized structures...\n')

    for name in low_e_names:
        update(name)

    print(f'Reading optimized structures...')

    coords, atomnos_list = [], []
    for name in low_e_names:
        mol = read_xyz(name.split('.')[0]+".xyz")
        atomnos_list.append(mol.atomnos)
        assert len(mol.atomnos) == len(atomnos_list[0]), f'Not all molecules have the same # of atoms! ({name})'
        coords.append(mol.atomcoords[0])

    coords = np.array(coords)
    graph = graphize(coords[0], atomnos_list[0])
    print(f'Removing similar structures...')
    coords = cl_similarity_refining(coords, atomnos_list[0], graph)

    with open(outname, "w") as f:
        for coord, atomnos in zip(coords, atomnos_list):
            write_xyz(coord, atomnos, f)

    print(f"Wrote {len(coords)} structures to {outname}.")