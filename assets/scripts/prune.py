import time
import sys

from firecode.utils import time_to_string, graphize, read_xyz, write_xyz
from firecode.pruning import prune_by_rmsd, prune_by_rmsd_rot_corr, prune_by_moment_of_inertia

def cl_similarity_refining(coords, atomnos, graph, moi=True, rmsd=True, rmsd_thr=0.25, verbose=False, payload=None):
    
        '''
        Removes structures that are too similar to each other (RMSD-based).
        payload: list of arrays that will be returned at the end, with each sequential maks applied to each element of the list

        '''

        if verbose:
            print('--> Similarity Processing')

        before = len(coords)

        if moi:

            if len(coords) <= 5E5:

                ### Prune based on the moment of inertia

                before3 = len(coords)
                t_start = time.perf_counter()

                coords, mask = prune_by_moment_of_inertia(coords, atomnos)

                if payload is not None:
                    payload = [elem[mask] for elem in payload]

                if before3 > len(coords):
                    print(f'Discarded {int(len([b for b in mask if not b]))} candidates for MOI similarity ({len([b for b in mask if b])} left, {time_to_string(time.perf_counter()-t_start)})')

        if rmsd and len(coords) <= 1E5:

            before1 = len(coords)
            t_start = time.perf_counter()

            # coords, mask = prune_conformers_rmsd(coords, atomnos, max_rmsd=self.options.rmsd, verbose=verbose)
            coords, mask = prune_by_rmsd(coords, atomnos, max_rmsd=rmsd_thr)

            if payload is not None:
                payload = [elem[mask] for elem in payload]

            if before1 > len(coords):
                print(f'Discarded {int(len([b for b in mask if not b]))} candidates for RMSD similarity ({len([b for b in mask if b])} left, {time_to_string(time.perf_counter()-t_start)})')

            ### Second step: again but symmetry-corrected (unless we have too many structures)

            if len(coords) <= 500:

                before2 = len(coords)
                t_start = time.perf_counter()

                coords, mask = prune_by_rmsd_rot_corr(coords, atomnos, graph, max_rmsd=rmsd_thr, logfunction=(print if verbose else None))

                if payload is not None:
                    payload = [elem[mask] for elem in payload]

                if before2 > len(coords):
                    print(f'Discarded {int(len([b for b in mask if not b]))} candidates for symmetry-corrected RMSD similarity ({len([b for b in mask if b])} left, {time_to_string(time.perf_counter()-t_start)})')


        if verbose and len(coords) == before:
            print(f'All structures passed the similarity check.{" "*15}')

        return coords, payload

if __name__ == '__main__':

    assert len(sys.argv) == 3, 'Usage: prune.py filename.xyz outname.xyz'
    _, filename, outname = sys.argv

    mol = read_xyz(filename)
    graph = graphize(mol.atomcoords[0], mol.atomnos)

    pruned, _ = cl_similarity_refining(mol.atomcoords, mol.atomnos, graph, verbose=True)

    with open(outname, "w") as f:
        for s in pruned:
            write_xyz(s, mol.atomnos, f)
