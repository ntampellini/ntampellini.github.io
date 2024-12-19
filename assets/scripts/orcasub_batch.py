import sys
import os
from InquirerPy import inquirer
from InquirerPy.base.control import Choice

partition = 'day'

def main(args):

    maxtime_string = {
        'day' : '0-24:00:00',
        'week': '7-00:00:00',

    }[partition]

    done = []

    if len(args) == 1:
        print(f"\n  Launch batches of ORCA jobs through orcasub.sh. Syntax:\n\n  " +
            "    python orcasub_batch.py conf*.xyz\n")
        quit()

    def get_procs(name):
        with open(name, "r") as f:
            lines = [line for line in f.readlines() if "nprocs" in line]
        if len(lines) != 1:
            return 16
        else:
            return int(lines[0].split()[1])

    def get_mem(name, safety_factor=1.333):
        with open(name, "r") as f:
            lines = [line for line in f.readlines() if "%maxcore" in line]
        if len(lines) != 1:
            return 8000
        else:
            return int(int(lines[0].split()[1])*safety_factor)

    for name in args[1:]:

        basename = name.split(".")[0]
        procs = get_procs(basename+".inp")
        mem = get_mem(basename+".inp")

        if basename not in done:

            if basename+".inp" in os.listdir() and basename+'.xyz' in os.listdir():

                os.system(f'sbatch -J 🐋⁶_ORCA_{os.path.basename(os.getcwd())}/{basename} --tasks-per-node {procs} --mem {int(procs*mem)} --partition {partition} -t {maxtime_string} orcasub_scratch.sh {basename}')
                print(f'Launched ORCA_{os.path.basename(os.getcwd())}/{basename} on {procs} cores / {float(procs*mem/1000):.2f} GB ({float(mem/1000):.2f} GB/core)')
                done.append(basename)

            else:

                f"Can't find {basename}.inp"

if __name__ == "__main__":

    partition = inquirer.select(
        message="Which partition would you like to run the jobs on?",
        choices=(
            Choice(value='day', name= 'day   - Max. duration 24 h'),
            Choice(value='week', name='week  - Max. duration 7 d'),
        ),
        default='day',
    ).execute()

    main(sys.argv)
