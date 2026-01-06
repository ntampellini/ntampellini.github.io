import os
import sys
from dataclasses import dataclass
from subprocess import getoutput

from rich.traceback import install
from utils import read_xyz, write_xyz
from InquirerPy import inquirer
from InquirerPy.base.control import Choice

install(show_locals=True, locals_max_length=None, locals_max_string=None, width=120)

def inquire_prev_to_last() -> bool:
    return inquirer.select(
                    message="What energy would you like to extract?",
                    choices=(
                        Choice(value=True, name='Previous to last energy (last at lower level)'),
                        Choice(value=False, name='Last energy (higher level)'),
                    ),
                    default=True,
                ).execute()


if len(sys.argv) == 1:
    print("\n  Updates ORCA input files by replacing them with the last step of the optimization trajectory. Syntax:\n\n" +
          "  python update.py conf*.xyz\n\n" + 
          "  conf*.xyz: base name of input geometry file(s)\n"
        )
    quit()

@dataclass
class Structure:

    basename: str
    prev_to_last: bool = False

    def __post_init__(self):
        self.charge, self.mult = getoutput(f'grep xyzfile {self.basename}.inp').split()[2:4]
        
        if self.prev_to_last:
            command = f'grep \"FINAL SINGLE POINT ENERGY\" {self.basename}.out | tail -2 | head -1'
        else:
            command = f'grep \"FINAL SINGLE POINT ENERGY\" {self.basename}.out | tail -1'
            
        self.el_energy = float(getoutput(command).rstrip(" Eh").split()[-1])
        self.freqs = self._read_freqs()

    def _read_freqs(self):
        n_atoms = int(getoutput(f'head {self.basename}.xyz -n 1'))
        n_freqs = 3 * n_atoms

        lines = getoutput(f'grep \" \+[0-9]\+: \+-*[0-9]\+\.[0-9]\+ cm\*\*-1\" {self.basename}.out').split('\n')
        freqs = [float(line.split()[1]) for line in lines]
        
        # if more than one set, only keep the last
        freqs = freqs[-n_freqs:]

        return freqs

    @property
    def data(self):
        return (f"charge: {self.charge}; multiplicity: {self.mult}; energy: {self.el_energy}; "
        + f"SMILES: xxx; frequencies: {self.freqs}; program: ORCA; DOI: xxx")

def update(name, prev_to_last=False):

    basename = name.split(".")[0]
    traj = basename + "_trj.xyz"
    xyzname = basename + ".xyz"
    outname = basename + ".out"
    files = os.listdir()

    if basename not in done:
        
        if traj in files:

            mol = read_xyz(traj)

            if outname in files:
                try:
                    title = Structure(basename).data
                except Exception as e:
                    print(e)
                    title = ""
                    pass
            else:
                title = ""

            with open(xyzname, "w") as f:
                write_xyz(mol.atomcoords[-1], mol.atomnos, f, title=title)

            print(f"Updated {xyzname}")
            
            done.append(name.split(".")[0])
            return True

        else:

            # update with the latest Compound job trajectory
            for i in range(4,0,-1):
                cmptraj = name.split(".")[0] + f"_Compound_{i}_trj.xyz"

                if cmptraj in os.listdir():

                    mol = read_xyz(cmptraj)

                    if outname in files:
                        try:
                            title = Structure(basename, prev_to_last=prev_to_last).data
                        except Exception as e:
                            print(e)
                            title = ""
                            pass
                    else:
                        title = ""

                    with open(xyzname, "w") as f:
                        write_xyz(mol.atomcoords[-1], mol.atomnos, f, title=title)

                    print(f"Updated {xyzname} from {cmptraj}")
                    
                    done.append(basename)
                    return True

            print(f"Can't find {traj} nor {cmptraj}")
            return False

if __name__ == "__main__":

    done = []
    inquire = True

    if inquire:
        prev_to_last = inquire_prev_to_last()
    else:
        prev_to_last = False

    for name in sys.argv[1:]:
        update(name, prev_to_last=prev_to_last)
        
    print(f"Updated {len(done)} structures.")