import sys
from compare import Job
from firecode.units import EH_TO_EV
import numpy as np
from firecode.utils import read_xyz
from subprocess import getoutput

def read_atomic_charges(job: Job, charge_type: str = "LOEWDIN") -> np.array:
    """Reads atomic charges of the desired type and returns them as an array."""

    lines = getoutput(f"grep -i \"{charge_type} ATOMIC\" {job.name} -A 500").split("\n")
    # print(f"ran command: grep -i \"{charge_type} ATOMIC\" {job.name} -A 500")
    # print(lines)
    output = []

    for line in lines:

        if line == "":
            break

        if ":" in line:
            output.append(float(line.split(":")[-1].split()[0]))

    print(f"{job.name}: read {len(output)} charges")
    return np.array(output)

def read_spin_densities(job: Job, spin_type: str = "LOEWDIN") -> np.array:
    """Reads spin densities of the desired type and returns them as an array."""

    lines = getoutput(f"grep -i \"{spin_type} ATOMIC\" {job.name} -A 500").split("\n")
    output = []

    for line in lines:

        if line == "":
            break

        if ":" in line:
            output.append(float(line.split(":")[-1].split()[1]))

    return np.array(output)


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("usage: python philicity.py Rdot.out R+.out R-.out")
        sys.exit(1)

    # radical = Job(sys.argv[1])
    # cation = Job(sys.argv[2])
    # anion = Job(sys.argv[3])

    anion, radical, cation = sorted([Job(f) for f in sys.argv[1:4]], key=lambda x: x.charge)

    assert radical.charge + 1 == cation.charge
    assert radical.charge - 1 == anion.charge

    print(f"Anion   : {anion.name}")
    print(f"Radical : {radical.name}")
    print(f"Cation  : {cation.name}")
    print()

    anion.read_energy()
    radical.read_energy()
    cation.read_energy()

    I = cation.get_comparison_energy() - radical.get_comparison_energy()
    A = anion.get_comparison_energy() - radical.get_comparison_energy()
    nu_hardness = I + A
    s_softness = 1 / nu_hardness

    anion_charges = read_atomic_charges(anion)
    radical_charges = read_atomic_charges(radical)
    cation_charges = read_atomic_charges(cation)

    fukui_plus = radical_charges - anion_charges
    fukui_minus = cation_charges - radical_charges
    fukui_zero = (fukui_plus + fukui_minus) / 2
    radical_philicity_indices = fukui_zero * s_softness

    omega_philicity = ((((I+A)**2)/(I-A))/8) * EH_TO_EV

    radical_spin_densities = read_spin_densities(radical)

    # print(f"{cation.name} --(+e^-)-> {radical.name} --(+e^-)-> {anion.name}")
    print()
    print(f"ω (global radical philicity)= {omega_philicity:.3f}")
    print(f"η (global hardness) = {nu_hardness:.3f}")
    print(f"S (global softness) = {s_softness:.3f}")
    print()

    sorting_ids = np.argsort(radical_philicity_indices)[::-1]
    radical.last_atoms = read_xyz(radical.name[:-4]+".xyz").atoms
    assert len(radical.atoms) == len(radical_charges), f"{len(radical.last_atoms)=} != {len(radical_charges)=}"

    print(f"Index  Symbol  P_k^0  Spin dens.")
    print("--------------------------------")
    for i, symbol, phil, dens in zip(
        sorting_ids,
        radical.last_atoms[sorting_ids],
        radical_philicity_indices[sorting_ids],
        radical_spin_densities[sorting_ids]):
        print(f"{i:3}  {symbol:2s}  {phil:+.4f}  {dens:+.4f}")





