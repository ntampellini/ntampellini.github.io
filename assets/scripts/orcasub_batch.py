import sys
import os

done = []
procs = 16

if len(sys.argv) == 1:
    print(f"\n  Launch batches of ORCA jobs through orcasub.sh. Syntax:\n\n  " +
          "python orcasub_batch.py conf*.xyz [n]\n\n  " +
          "n: optional, overrides the default value of 16 for cores\n")
    quit()

try:
    if int(sys.argv[-1]) in (1, 2, 4, 8, 16, 32, 64):
        *sys.argv, procs = sys.argv
        
except ValueError:
    pass

print(f'--> Running jobs on {procs} cores each')

for name in sys.argv[1:]:

    basename = name.split(".")[0]

    if basename not in done:
        
        if basename+".inp" in os.listdir() and basename+'.xyz' in os.listdir():

            os.system(f'sbatch orcasub.sh {basename} {procs}')
            print(f'Launched {basename}')            
            done.append(basename)

        else:

            f"Can't find {basename}.inp"
