import sys
import os

done = []

if len(sys.argv) == 1:
    print(f"\n  Launch batches of ORCA jobs through orcasub.sh. Syntax:\n\n  " +
          "python orcasub_batch.py conf*.xyz\n\n")
    quit()

for name in sys.argv[1:]:

    basename = name.split(".")[0]

    if basename not in done:

        if basename+".inp" in os.listdir() and basename+'.xyz' in os.listdir():

            os.system(f'sbatch -J ORCA_{os.path.basename(os.getcwd())}/{basename} orcasub_scratch.sh {basename}')
            print(f'Launched ORCA_{os.path.basename(os.getcwd())}/{basename}')
            done.append(basename)

        else:

            f"Can't find {basename}.inp"
