import os
import sys
from subprocess import getoutput

look_for = "Final Gibbs free energy"
completed, incomplete = [], []


for name in sys.argv[1:]:
    print(f'Reading {name}...', end="\r")
    try:
        energy = float(getoutput(f'grep \"{look_for}\" {name} | tail -1').rstrip(" Eh").split()[-1])
        completed.append(name)
    
    except IndexError:
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
        main([None, *incomplete])