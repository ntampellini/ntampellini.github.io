import os
import sys

names = []
home = os.getcwd()

for folder in sys.argv[1:]:
    os.chdir(os.path.join(home, folder))
             
    for filename in os.listdir():
        if filename.endswith(".xyz"):
            # names.append(filename)
            with open(filename, "r") as f:
                lines = f.readlines()
                n_structs = lines.count(lines[0])

            print(f'{filename:20} {n_structs:3} structures')

# print(f'Found {len(names)} .xyz files.')

# for name in names:
    # mol = ccread(name)