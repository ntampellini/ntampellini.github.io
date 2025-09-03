import sys
import os
from InquirerPy import inquirer
from InquirerPy.base.control import Choice

def get_procs(name, default=1):
    with open(name, "r") as f:
        lines = [line for line in f.readlines() if "nprocs" in line]
    if len(lines) != 1:
        return default
    else:
        return int(lines[0].split()[1])

def get_mem(name, safety_factor=None, default=8000):
    
    if safety_factor is None:
        safety_factor = 1.333
    
    with open(name, "r") as f:
        lines = [line for line in f.readlines() if "%maxcore" in line]
    if len(lines) != 1:
        return default
    else:
        return int(int(lines[0].split()[1])*safety_factor)

def main(rootnames, priority=False):

    if len(rootnames) == 0:
        print(f"\n  Launch batches of ORCA jobs through orcasub.sh. Syntax:\n\n  " +
            "    python orcasub_batch.py conf*.xyz\n")
        quit()

    priority_string = '--priority' if priority else ''

    partition = inquirer.select(
        message="Which partition would you like to run the jobs on?",
        choices=(
            Choice(value='day', name= 'day   - Max. duration 24 h'),
            Choice(value='week', name='week  - Max. duration 7 d'),
        ),
        default='day',
    ).execute()

    safety_factor = inquirer.text(
        message="Memory safety factor? (MEM requested/MEM allowed ORCA to use)",
        default="1.333",
        filter=float,
        validate=lambda x: float(x) > 1.25,
        invalid_message='Minimum safety factor is 1.25'
    ).execute()

    maxtime_string = {
        'day' : '0-24:00:00',
        'week': '7-00:00:00',
    }[partition]

    done = []
    for name in rootnames:

        basename = name.split(".")[0]
        procs = get_procs(basename+".inp")
        mem = get_mem(basename+".inp", safety_factor=safety_factor)

        if basename not in done:

            if basename+".inp" in os.listdir():
                # and basename+'.xyz' in os.listdir():

                os.system((f'sbatch -J 🐋⁶_ORCA_{os.path.basename(os.getcwd())}/{basename} --tasks-per-node {procs} ' + 
                          f'--mem {int(procs*mem)} --partition {partition} -t {maxtime_string} orcasub_scratch.sh {basename}{priority_string}'))
                print(f'Launched ORCA_{os.path.basename(os.getcwd())}/{basename} on {procs} cores / {float(procs*mem/1000):.2f} GB ({float(mem/1000):.2f} GB/core) on {partition}')
                done.append(basename)

            else:
                f"Can't find {basename}.inp"

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='description')
    parser.add_argument('filenames', action='store', nargs="*", help='filenames (rootname.*)', default=None)
    parser.add_argument('-p', '--priority', action="store_true", help='priority keyword', required=False, default=False)

    args = parser.parse_args(sys.argv[1:])

    if args.priority:
        print('--> Priority flag: running jobs in priority mode')

    main(args.filenames, args.priority)