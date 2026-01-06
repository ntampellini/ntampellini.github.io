####################################################
options = {
    "method" : "GFN2-XTB",
    "solvent" : "acetonitrile",
    # "charge" : 0,
}
####################################################

from optparse import OptionParser
from firecode.calculators._xtb import xtb_opt, xtb_get_free_energy
from firecode.utils import read_xyz
import os 
import sys

for option, value in options.items():
    print(f"--> {option} = {value}")

os.chdir(os.getcwd())
names, final_names = [], []

if len(sys.argv) < 2:
    for f in os.listdir():
        if os.path.isfile(f) and f.split('.')[1] == 'xyz':
            names.append(f)
else:
    for name in sys.argv[1:]:
        names.append(name)

energies, free_energies = [], []

for i, name in enumerate(names):

    charge = name.count('+') - name.count('-')

    try:
        data = read_xyz(name)
            
        print(f'Calculating SP energy on {name} - {i+1} of {len(names)} ({options["method"]})')
        coords, energy, success = xtb_opt(data.atomcoords[0], data.atomnos,
            method=options["method"], solvent=options["solvent"], charge=charge, opt=False)
        energies.append(energy)

        print(f'Calculating Free Energy on {name} - {i+1} of {len(names)} ({options["method"]})')
        free_energy = xtb_get_free_energy(data.atomcoords[0], data.atomnos, method=options["method"],
            solvent=options["solvent"], charge=charge, sph=True, debug=True)
        free_energies.append(free_energy)

        final_names.append(name)

    except IndexError:
        pass

print("> NAME                    EE(GFN2)           G-E(el)       G(GFN2)       Rel. G (kcal/mol)")
print("-"*100)

min_fe = min(free_energies)

for name, energy, free_energy in zip(final_names, energies, free_energies):
# for name, energy, free_energy in zip(final_names, energies, energies):
    print("> {:20}    {:10}    {:10}    {:10}    {:6}".format(name,
                                                                round((energy) / 627.5096080305927, 6),
                                                                round((free_energy-energy) / 627.5096080305927, 6),
                                                                round((free_energy) / 627.5096080305927, 6),
                                                                round((free_energy-min_fe), 2)))
