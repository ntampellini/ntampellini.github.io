import sys
import os

def main(args):

    done = []

    if len(args) == 1:
        print(f"\n  Launch batches of ORCA jobs through orcasub.sh. Syntax:\n\n  " +
            "    python orcasub_batch.py conf*.xyz [32]\n\n" +
            "  32: if specified, overrides the default value for cores (16)")
        quit()

    def get_procs(name):
        with open(name, "r") as f:
            lines = [line for line in f.readlines() if "nprocs" in line]
        if len(lines) != 1:
            return 16
        else:
            return int(lines[0].split()[1])

    for name in args[1:]:

        basename = name.split(".")[0]
        procs = get_procs(basename+".inp")

        if basename not in done:

            if basename+".inp" in os.listdir() and basename+'.xyz' in os.listdir():

                os.system(f'sbatch -J ORCA_{os.path.basename(os.getcwd())}/{basename} --tasks-per-node {procs} --mem {int(procs*8)}G orcasub_scratch.sh {basename}')
                print(f'Launched ORCA_{os.path.basename(os.getcwd())}/{basename} on {procs} cores / {int(procs*8)} GB')
                done.append(basename)

            else:

                f"Can't find {basename}.inp"

if __name__ == "__main__":
    main(sys.argv)