import argparse
from subprocess import getoutput

import numpy as np
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from rich.traceback import install

from compare import Job, assert_homogeneous_temps
from utils import tcolors

EH_TO_KCAL = 627.5096080305927
KCAL_TO_EV = 0.0433641153087705
R = 0.001985877534

install(show_locals=True)

class Reaction:
    def __init__(self, reactants, products):
        self.reactants: list[Job] = reactants
        self.products: list[Job] = products

        r_atoms = sum([job.natoms for job in self.reactants])
        p_atoms = sum([job.natoms for job in self.products])
        if r_atoms == p_atoms:
            print(f'\n--> Number of atoms check passed: {r_atoms} atoms in both reactants and products.')
        else:
            raise ValueError(f"Number of atoms in reactants ({r_atoms}) does not match number of atoms in products ({p_atoms}).")
        
        r_charge = sum([job.charge for job in self.reactants])
        p_charge = sum([job.charge for job in self.products])
        if r_charge == p_charge:
            print(f'--> Total charge check passed: {r_charge} in both reactants and products.')
        else:
            raise ValueError(f"Total charge in reactants ({r_charge}) does not match total charge in products ({p_charge}).")

        # Now check vibrations and print warnings
        vibs_ok = True
        for reactant in self.reactants:
            if not self.check_vibs(reactant):
                vibs_ok = False

        for product in self.products:
            if not self.check_vibs(product):
                vibs_ok = False

        if vibs_ok:
            print(f"--> Number of negative vibrations in all species matches expectations (0 unless name starts with \"ts\").")

    def check_vibs(self, job: Job) -> bool:
        """Check vibrations and warn if they look suspect."""

        n_freqs = 3 * job.natoms
        is_ts = (getoutput(f'grep -i \"\\!.*optts\" {job.freqfile}') != "") or job.name.lower().startswith("ts")

        lines = getoutput(f'grep \" \\+[0-9]\\+: \\+-*[0-9]\\+\\.[0-9]\\+ cm\\*\\*-1\" {job.freqfile}').split('\n')
        freqs = [float(line.split()[1]) for line in lines]

        if len(freqs) < n_freqs:
            raise Exception(f"No frequencies read for {job.freqfile}")
        
        freqs = freqs[-n_freqs:]
        neg_freqs = len([f for f in freqs if f < 0])
        
        if is_ts and neg_freqs == 1:
            return True
        
        if not is_ts and neg_freqs == 0:
            return True
        
        print(f"{tcolors.BOLD}{tcolors.YELLOW}--> WARNING: {job.freqfile} has {neg_freqs} "
              f"negative frequenc{"y" if neg_freqs == 1 else "ies"}. Make sure this is "
              f"the expected number.{tcolors.ENDC}")
        
        return False


    @property
    def delta_free_energy(self):
        if not self.reactants or not self.products:
            raise ValueError("Both reactants and products must be defined to calculate delta energy.")
        reactant_energy = sum([job.free_energy for job in self.reactants])
        product_energy = sum([job.free_energy for job in self.products])
        return product_energy - reactant_energy
    
    @property
    def delta_enthalpy(self):
        if not self.reactants or not self.products:
            raise ValueError("Both reactants and products must be defined to calculate delta energy.")
        reactant_energy = sum([job.enthalpy for job in self.reactants])
        product_energy = sum([job.enthalpy for job in self.products])
        return product_energy - reactant_energy
        
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("outfiles", help="ORCA Outfiles to include in the interactive menu.", action="store", nargs="*")
    parser.add_argument("-r", "--reactants", help="ORCA Outfiles to use as reactants.", action="store", nargs="*", required=False)
    parser.add_argument("-p", "--products", help="ORCA Outfiles to use as products.", action="store", nargs="*", required=False)
    args = parser.parse_args()

    if args.reactants and args.products:
        reactants = args.reactants
        products = args.products

    else:
        outnames = getoutput("ls *.out").split('\n')
        outnames += args.outfiles

        outnames = sorted(list(set(outnames))) # Removes duplicates and sorts
        choices = [Choice(value=f, name=f) for f in outnames if f]

        reactants = inquirer.checkbox(
            choices=choices,
            message="Select reactant files (spacebar to toggle, enter to confirm):",
            disabled_symbol='⬡',
            enabled_symbol='⬢',
        ).execute()

        products = inquirer.checkbox(
            choices=choices,
            message="Select reactant files (spacebar to toggle, enter to confirm):",
            disabled_symbol='⬡',
            enabled_symbol='⬢',
        ).execute()

    avail_gcorrs = []
    if "GIBBS" in getoutput(f'grep GIBBS {reactants[0]}'):
        avail_gcorrs.append(Choice(value=".", name=f'./    This folder   - read free energy from this folder.'))

    if "GIBBS" in getoutput(f'grep GIBBS ../{reactants[0]}'):
        avail_gcorrs.append(Choice(value="..", name=f'../   Parent folder - read free energy from the parent folder.'))

    avail_gcorrs.append(Choice(value=None,    name=f'?     Other folder  - choose another folder to read G(corr) values.'))

    freqdir = inquirer.select(
        message=f"Which folder would you like to extract G/H(corr) values from?",
        choices=avail_gcorrs,
        default=avail_gcorrs[0].value,
    ).execute()

    reactant_jobs = []
    for filename in reactants:
        job = Job(filename, energy_mode="G", freqdir=freqdir)
        job.read_electronic_energy()
        job.read_enthalpy()
        job.read_free_energy()
        reactant_jobs.append(job)

    assert_homogeneous_temps(reactant_jobs)

    product_jobs = []
    for filename in products:
        job = Job(filename, energy_mode="G", freqdir=freqdir)
        job.read_electronic_energy()
        job.read_enthalpy()
        job.read_free_energy()
        product_jobs.append(job)

    T = assert_homogeneous_temps(product_jobs)

    reaction = Reaction(reactant_jobs, product_jobs)
    print(f'--> Temperature check passed: {T:.2f} K ({T-273.15:+.1f} °C) in all reactants and products.')

    dG = reaction.delta_free_energy * EH_TO_KCAL
    dG_eV = dG * KCAL_TO_EV
    dH = reaction.delta_enthalpy * EH_TO_KCAL
    minus_TdS = dG - dH
    dS = -minus_TdS / T
    K_eq = np.exp(-dG/(R*T))
    K_eq_inv = np.exp(dG/(R*T))

    print(f'\n{rxn_color(-dG)}{" + ".join(reactants)}{tcolors.ENDC}{arrow(dG)}{rxn_color(dG)}{" + ".join(products)}{tcolors.ENDC}\n')
    print(f'Reaction Free Energy ({T:.2f} K):    ΔG° = {value_color(dG)}{dG:+.2f}{tcolors.ENDC} kcal/mol ({dG_eV:+.2f} eV)')
    print(f'Reaction Enthalpy    ({T:.2f} K):    ΔH° = {value_color(dH)}{dH:+.2f}{tcolors.ENDC} kcal/mol')
    print(f'Reaction -TΔS        ({T:.2f} K):  -TΔS° = {value_color(minus_TdS)}{minus_TdS:+.2f}{tcolors.ENDC} kcal/mol')
    print(f'Reaction Entropy     ({T:.2f} K):    ΔS° = {value_color(dS)}{dS*1000:+.1f}{tcolors.ENDC} cal/(mol*K) [e.u.]\n')
    print(f'K_eq                 ({T:.2f} K):   K_eq = {K_eq:.3e}')
    print(f'log10(K_eq)          ({T:.2f} K): log(K) = {np.log10(K_eq):.3f}')
    print(f'(K_eq)^-1            ({T:.2f} K):  K_inv = {K_eq_inv:.3e}')

def rxn_color(value: float, thr: float = 1.0) -> str:
    """Returns the appropriate tcolor value."""
    if value > thr:
        return tcolors.ENDC
    
    if value < -thr:
        return tcolors.BOLD
    
    return tcolors.ENDC

def value_color(value: float, thr: float = 1.0) -> str:
    """Returns the appropriate tcolor value."""
    if value > thr:
        return tcolors.BOLD + tcolors.RED
    
    if value < -thr:
        return tcolors.BOLD + tcolors.GREEN
    
    return tcolors.ENDC

def arrow(dG: float, thr: float = 1.5) -> str:
    if dG > thr:
        return " <-- "

    if dG < -thr:
        return " --> "
    
    return " <=> "


if __name__ == '__main__':
    main()