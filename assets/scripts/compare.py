from subprocess import getoutput
import sys

if len(sys.argv) == 1:
    print(f"\n  Compare completed ORCA jobs electronic/free energy of groups of jobs and prints " +
           "collected energy values. Syntax:\n\n" +
           "  python compare.py conf*.out [g]\n\n" + 
           "  conf*.out: base name with generic asterisk descriptor (group of output files)\n" +
           "  g: optional, uses free energy for the comparison (needs freq calculations).\n")
    quit()

g = False

if "g" in sys.argv:
    g = True
    look_for = "Final Gibbs free energy"
    sys.argv.remove("g")
else:
    look_for = "FINAL SINGLE POINT ENERGY"

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

for name, energy in zip(names, energies):
    converged = True if getoutput(f'grep HURRAY {name}') != '' else False
    iterations = "conv" if converged else getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {name} -c')
    print(f'{name} - {round((energy-min_e)*627.509608030593, 3)} - [{iterations}{"" if converged else " iterations"}]')

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