import sys
from subprocess import getoutput

from compare import Job, assert_homogeneous_temps
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from rich.traceback import install

EH_TO_KCAL = 627.5096080305927

install(show_locals=True)

class Reaction:
    def __init__(self, reactants, products):
        self.reactants: list[Job] = reactants
        self.products: list[Job] = products

        r_atoms = sum([job.natoms for job in self.reactants])
        p_atoms = sum([job.natoms for job in self.products])
        if r_atoms == p_atoms:
            print(f'--> Number of atoms check passed: {r_atoms} atoms in both reactants and products.')
        else:
            raise ValueError(f"Number of atoms in reactants ({r_atoms}) does not match number of atoms in products ({p_atoms}).")
        
        r_charge = sum([job.charge for job in self.reactants])
        p_charge = sum([job.charge for job in self.products])
        if r_charge == p_charge:
            print(f'--> Total charge check passed: {r_charge} in both reactants and products.')
        else:
            raise ValueError(f"Total charge in reactants ({r_charge}) does not match total charge in products ({p_charge}).")

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
    
    outnames = getoutput("ls *.out").split('\n')
    outnames += sys.argv[1:]
    outnames = sorted(list(set(outnames))) # Remove duplicates adn sorts
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

    print(f'\n{" + ".join(reactants)} --> {" + ".join(products)}\n')
    print(f'Reaction Free Energy ({T:.2f} K): ΔG° = {reaction.delta_free_energy*EH_TO_KCAL:.2f} kcal/mol')
    print(f'Reaction Enthalpy    ({T:.2f} K): ΔH° = {reaction.delta_enthalpy*EH_TO_KCAL:.2f} kcal/mol\n')
    
if __name__ == '__main__':
    main()