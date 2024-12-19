import os 

# peptide, sub
pairs = [
    (1, 6),
    (3, 6),

    (2, 9),
    (4, 9),

    (3, 4),
    (4, 4),

    (4, 10),
    (5, 10),

    (4, 11),
    # (6, 11),

    (5, 12),
    (7, 12),
]

for peptide, substrate in pairs:
    if os.path.isdir(f'./sub{substrate}'):

        if substrate == 12:
            os.system(f'python /gpfs/gibbs/project/miller/nt383/scripts/compare.py sub{substrate}/p{peptide}/*/ent/gs_*/gs*.out g x=/gpfs/gibbs/project/miller/shared/nt383/diols/completed/{peptide}e{substrate}.xyz')
            os.system(f'python /gpfs/gibbs/project/miller/nt383/scripts/compare.py sub{substrate}/p{peptide}/*/report/gs_*/gs*.out g x=/gpfs/gibbs/project/miller/shared/nt383/diols/completed/{peptide}r{substrate}.xyz')

        else:
            os.system(f'python /gpfs/gibbs/project/miller/nt383/scripts/compare.py sub{substrate}/p{peptide}/*/ent/gs/gs*.out g x=/gpfs/gibbs/project/miller/shared/nt383/diols/completed/{peptide}e{substrate}.xyz')
            os.system(f'python /gpfs/gibbs/project/miller/nt383/scripts/compare.py sub{substrate}/p{peptide}/*/report/gs/gs*.out g x=/gpfs/gibbs/project/miller/shared/nt383/diols/completed/{peptide}r{substrate}.xyz')


    else:
        print(f"Could not find folder \"sub{substrate}\"")
