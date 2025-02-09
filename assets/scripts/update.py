from utils import read_xyz, write_xyz
import sys
import os

done = []

if len(sys.argv) == 1:
    print("\n  Updates ORCA input files by replacing them with the last step of the optimization trajectory. Syntax:\n\n" +
          "  python update.py conf*.xyz\n\n" + 
          "  conf*.xyz: base name of input geometry file(s)\n"
        )
    quit()

def update(name):

    traj = name.split(".")[0] + "_trj.xyz"
    inp = name.split(".")[0] + ".xyz"

    if name.split(".")[0] not in done:
        
        if traj in os.listdir():

            mol = read_xyz(traj)

            with open(inp, "w") as f:
                write_xyz(mol.atomcoords[-1], mol.atomnos, f)

            print(f"Updated {inp}")
            
            done.append(name.split(".")[0])
            return True

        else:

            # update with the latest Compound job trajectory
            for i in range(4,0,-1):
                cmptraj = name.split(".")[0] + f"_Compound_{i}_trj.xyz"

                if cmptraj in os.listdir():

                    mol = read_xyz(cmptraj)

                    with open(inp, "w") as f:
                        write_xyz(mol.atomcoords[-1], mol.atomnos, f)

                    print(f"Updated {inp} from {cmptraj}")
                    
                    done.append(name.split(".")[0])
                    return True

            print(f"Can't find {traj} nor {cmptraj}")
            return False

if __name__ == "__main__":

    for name in sys.argv[1:]:
        update(name)
        
    print(f"Updated {len(done)} structures.")