import os
import sys
from subprocess import getoutput

# look_for = "Final Gibbs free energy"
look_for = "ORCA TERMINATED NORMALLY"
completed, incomplete = [], []


for name in sys.argv[1:]:

    basename = name.split(".")[0]
    outname = basename + ".out"

    print(f'Reading {basename}...', end="\r")
    
    if outname in os.listdir() and getoutput(f'grep \"{look_for}\" {outname}') != "":
        completed.append(name)

    else:
        incomplete.append(name)

if len(incomplete) == 0:
    print(f"All {len(sys.argv[1:])} output files look completed (contain string \"{look_for}\" at least one time).")

else:
    print(f"Found {len(incomplete)}/{len(sys.argv[1:])} incomplete jobs:")
    for name in incomplete:
        basename = name.split(".")[0]
        print(f"  - {basename}")

    print()
    
    answer = None
    while answer not in ("y", "n"):
        answer = input(f"Do you want to resubmit these {len(incomplete)} jobs? [y]/n:")
        answer = "y" if answer == "" else answer

    if answer == "y":

        from orcasub_batch import main
        main(incomplete)