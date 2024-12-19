import sys

import numpy as np
from firecode.algebra import all_dists
from firecode.pt import pt
from firecode.rmsd import rmsd_and_max_numba
from firecode.torsion_module import (_get_hydrogen_bonds, _get_torsions,
                                     _is_nondummy, get_double_bonds_indices,
                                     rotationally_corrected_rmsd_and_max)
from firecode.utils import graphize, read_xyz
from networkx import connected_components


def rmsd_rot_corr_easy(ref, tgt, graph, atomnos, logfunction=None):

    # center structures
    ref -= ref.mean(axis=0)
    tgt -= tgt.mean(axis=0)

    # get the number of molecular fragments
    subgraphs = list(connected_components(graph))

    # if they are more than two, give up on pruning by rot corr rmsd
    if len(subgraphs) > 2:
        raise NotImplementedError
    
    # if they are two, we can add a fictitious bond between the closest
    # atoms on the two molecular fragment in the provided graph, and
    # then removing it before returning
    if len(subgraphs) == 2:
        subgraphs =  [list(set) for set in connected_components(graph)]
        all_dists_array = all_dists(ref[list(subgraphs[0])], ref[list(subgraphs[1])])
        min_d = np.min(all_dists_array)
        s1, s2 = np.where(all_dists_array == min_d)
        i1, i2 = subgraphs[0][s1[0]], subgraphs[1][s2[0]]
        graph.add_edge(i1, i2)
        
    # add hydrogen bonds to molecular graph 
    hydrogen_bonds = _get_hydrogen_bonds(ref, atomnos, graph)
    for hb in hydrogen_bonds:
        graph.add_edge(*hb)

    # get all rotable bonds in the molecule, including dummy rotations
    torsions = _get_torsions(graph,
                            hydrogen_bonds=_get_hydrogen_bonds(ref, atomnos, graph),
                            double_bonds=get_double_bonds_indices(ref, atomnos),
                            keepdummy=True,
                            mode='symmetry')

    # only keep dummy rotations (checking both directions)
    torsions = [t for t in torsions if not (
                                _is_nondummy(t.i2, t.i3, graph) and (
                                _is_nondummy(t.i3, t.i2, graph)))]

    # since we only compute RMSD based on heavy atoms, discard quadruplets that involve hydrogen atoms
    torsions = [t for t in torsions if 1 not in [atomnos[i] for i in t.torsion]]

    # get torsions angles
    angles = [t.get_angles() for t in torsions]

    # Used specific directionality of torsions so that we always rotate the dummy portion (the one attached to the last index)
    torsions = [list(t.torsion) if _is_nondummy(t.i2, t.i3, graph) else list(reversed(t.torsion)) for t in torsions]
  

    # Print out torsion information
    if logfunction is not None:
        logfunction('\n >> Dihedrals considered for subsymmetry corrections:')
        for i, (torsion, angle) in enumerate(zip(torsions, angles)):
            logfunction(' {:2s} - {:21s} : {}{}{}{} : {}-fold'.format(
                                                                str(i+1),
                                                                str(torsion),
                                                                pt[atomnos[torsion[0]]].symbol,
                                                                pt[atomnos[torsion[1]]].symbol,
                                                                pt[atomnos[torsion[2]]].symbol,
                                                                pt[atomnos[torsion[3]]].symbol,
                                                                len(angle)))
        logfunction("\n")

    result = rotationally_corrected_rmsd_and_max(ref, tgt, atomnos, torsions=torsions, graph=graph, angles=angles)

    # remove the extra bond in the molecular graph
    if len(subgraphs) == 2:
        graph.remove_edge(i1, i2)

    return result

if __name__ == '__main__':
    mol1 = read_xyz(sys.argv[1])
    mol2 = read_xyz(sys.argv[2])
    rmsd, maxdev = rmsd_and_max_numba(mol1.atomcoords[0], mol2.atomcoords[0], center=True)
    print()
    print(f'RMSD    : {rmsd:.2f} A')
    print(f'Max Dev.: {maxdev:.2f} A')

    graph = graphize(mol1.atomcoords[0], mol1.atomnos)
    rmsd, maxdev = rmsd_rot_corr_easy(mol1.atomcoords[0], mol2.atomcoords[0], graph=graph, atomnos=mol1.atomnos)
    print()
    print(f'(rot. corr.) RMSD    : {rmsd:.2f} A')
    print(f'(rot. corr.) Max Dev.: {maxdev:.2f} A')