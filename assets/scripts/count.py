import os

from cclib.io import ccread

import sys



names = []

home = os.getcwd()
total = 0



for folder in sys.argv[1:]:

    os.chdir(os.path.join(home, folder))

             

    for filename in os.listdir():

        if filename.endswith(".xyz"):

            # names.append(filename)

            with open(filename, "r") as f:

                lines = f.readlines()

                n_structs = lines.count(lines[0])
                total += n_structs



            print(f'{filename:20} {n_structs:3} structures')


print(f'Total: {total} structures')
# print(f'Found {len(names)} .xyz files.')



# for name in names:

    # mol = ccread(name)

    
