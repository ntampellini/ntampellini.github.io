import os
import sys
from subprocess import getoutput

from InquirerPy import inquirer
from prettytable import PrettyTable

from utils import tcolors

# look_for = "Final Gibbs free energy"
look_for = "ORCA TERMINATED NORMALLY"
completed, incomplete = [], []

assert len(sys.argv) > 1, "Please provide at least one file as argument."

running = [line.split("/")[-1] for line in getoutput("squeue --me --format \"%100j\"").split()[1:]]
table = PrettyTable(["Job basename", "Status"])

for name in sys.argv[1:]:

    basename = name.split(".")[0]
    outname = basename + ".out"
    
    if outname in os.listdir() and getoutput(f'grep \"{look_for}\" {outname}') != "":
        completed.append(name)
        table.add_row([basename, f"{tcolors.BOLD+tcolors.GREEN}COMPLETED{tcolors.ENDC}"])

    elif basename not in running:
        incomplete.append(name)
        table.add_row([basename, f"{tcolors.BOLD+tcolors.RED}FAILED{tcolors.ENDC}"])

    else:
        table.add_row([basename, f"{tcolors.BOLD+tcolors.YELLOW}RUNNING{tcolors.ENDC}"])

print(table.get_string())

if len(incomplete) == 0:
    print(f"\nAll {len(sys.argv[1:])} jobs look completed (output contains \"{look_for}\") or are running.")

else:
    print(f"\nFound {len(incomplete)}/{len(sys.argv[1:])} incomplete jobs not currently running:")
    for name in incomplete:
        basename = name.split(".")[0]
        print(f"  - {basename}")

    print()
    
    if inquirer.confirm(
        message=f"Do you want to resubmit these {len(incomplete)} jobs?",
        default=True,
    ).execute():

        from orcasub_batch import main
        main(incomplete)