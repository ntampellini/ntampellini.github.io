import sys
from cclib.io import ccread
import numpy as np

print("--> CREST Constraint file writer <--\n")

if len(sys.argv) != 2:
    print("Specify one input file as an argument")

constraints, constrained = [], []
while True:

    data = input("Constrained indices: ").split()

    if data == []:
        break

    try:
        indices = [int(x) for x in data]
        assert len(indices) > 1
        constraints.append(indices)
        for i in indices:
            constrained.append(i)
    except Exception:
        pass

n_atoms = len(ccread(sys.argv[1]).atomnos)
types = {
    2 : "distance",
    3 : "angle",
    4 : "dihedral",
}

print(f"Specified {len(constraints)} constraints, total of {n_atoms} atoms")
    
with open("constraints.inp", "w") as f:
    f.write("$constrain\n  force constant=0.25\n")
    for indices in constraints:
        c_type = types[len(indices)]
        aug_ids = str([i+1 for i in indices]).strip("[").strip("]")
        f.write(f"  {c_type}: {aug_ids}, auto\n")

    f.write("$metadyn\n  atoms: ")

    # write atoms that need to be moved during metadynamics (all but constrained)
    l = np.array([i for i in range(1,n_atoms+1) if i-1 not in constrained])
    s = ""
    while len(l) > 2:
        i = next((i for i, _ in enumerate(l[:-2]) if l[i+1]-l[i]>1), len(l)-1)
        if l[0] == l[i]:
            s += f"{l[0]},"
        else:
            s += f"{l[0]}-{l[i]},"
        l = l[i+1:]

    # remove final comma
    f.write(s[:-1])

    f.write("\n$end")

print(f"Wrote constraints.inp file.")