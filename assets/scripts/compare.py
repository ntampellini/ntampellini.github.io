import os
import sys
from subprocess import getoutput

from firecode.utils import graphize, read_xyz, write_xyz, suppress_stdout_stderr

scratchdir = '/vast/palmer/scratch/miller/nt383'
EH_TO_KCAL = 627.509608030593

class Job:
    def __init__(self, filename, free_energy=False):
        self.name = filename

        if free_energy:
            self.compare_via = 'free_energy'
        else:
            self.compare_via = 'electronic_energy'

    def __gt__(self, other):
        return getattr(self, self.compare_via) > getattr(other, other.compare_via)
    
    def get_comparison_energy(self):
        return getattr(self, self.compare_via)

def compare(argv):
    '''
    '''

    if len(argv) == 1:
        print("\n  Compare completed ORCA jobs electronic/free energy of groups of jobs and prints " +
            "collected energy values. Syntax:\n\n" +
            "  python compare.py conf*.out [g] [x=folder/[file.xyz]]\n\n" + 
            "  conf*.out: base name with generic asterisk descriptor (group of output files)\n" +
            "  g: optional, uses free energy for the comparison (needs freq calculations).\n" +
            "  x: extract conformers from the specified files, cutting at [thr] kcal/mol and removing duplicates.\n" +
            "     Writes structures to folder/file.xyz\n" +
            "  thr: energy threshold for structure extraction, in kcal/mol [usage: thr=10]" +
            "  composite: builds Free energy using a single point energy calculation in the current\n" +
            "     folder and a frequency calculation in the [freq] folder, which is \"../\" by default.\n" +
            "  freq: specify the frequency folder for [composite] [usage: freq=../]")
        quit()

    ### If checking scratch, append the running job full names to argv

    if "scratch" in argv:
        argv.remove("scratch")
        
        cwd = os.getcwd()

        scratchnames = getoutput(f'ls {scratchdir}').split()
        
        pids = [line.split()[0] for line in getoutput("squeue --me").split('\n')[1:]]
        for pid in pids:
            try:
                workdir = [line.split('=')[1] for line in getoutput(f"scontrol show job {pid}").split("\n") if 'WorkDir' in line][0]
                if workdir == cwd:
                    jobdir = next(name for name in scratchnames if pid in name)
                    jobname = jobdir.split('-')[1].strip(pid).strip('_')
                    fullname = scratchdir+"/"+jobdir+'/'+jobname+'.out'
                    argv.append(fullname)

                    with suppress_stdout_stderr():
                        # copy traj from scratch
                        os.system(f'cp {scratchdir}/{jobdir}/{jobname}_trj.xyz .')
                        os.system(f'cp {scratchdir}/{jobdir}/{jobname}_Compound_*_trj.xyz .')
                        os.system(f'cp {scratchdir}/{jobdir}/{jobname}.property.txt .')
                        # os.system(f'cp {scratchdir}/{jobdir}/{jobname}.out .')

            except StopIteration:
                continue

    ### If interested in a composite method, see if the
    ### user specified where to extract gcorrs

    if "freq" in [kw.split("=")[0] for kw in argv]:
        freqdir = next((kw.split("=")[-1] for kw in argv if "freq" in kw))
        argv.remove(f"freq={freqdir}")
    else:
        freqdir = '..'

    g = False
    composite = False
    if "composite" in argv:
        argv.remove("composite")
        composite = True
        g = True
        folderstring = 'parent folder' if freqdir == '..' else freqdir
        print(f'--> COMPOSITE: extracting gcorrs from {folderstring}')

    ### If we are interested in free energy, set appropriate variables

    if "g" in argv:
        g = True
        argv.remove("g")

    extract_to_files = False
    energy_thr = 10.0

    ### Set the energy threshold for file extraction

    if "thr" in [kw.split("=")[0] for kw in argv]:
        energy_thr = next((kw.split("=")[-1] for kw in argv if "thr=" in kw))
        print(f'-> set thr to {energy_thr} kcal/mol')
        try:
            argv.remove(f"thr={energy_thr}")
        except ValueError:
            raise ValueError(f'Trying to delete non-existent\"thr={energy_thr}\"')
        
        energy_thr = float(energy_thr)

    # Specify where to extract files - if a filename is provided, structures
    # will be extracted to a single file, otherwise the same filenames of
    # initial structures will be used

    if "x" in [kw.split("=")[0] for kw in argv]:
        outname = next((kw.split("=")[-1] for kw in argv if "=" in kw))
        try:
            argv.remove(f"x={outname}")
        except ValueError:
            raise ValueError(f'Trying to delete non-existent\"x={outname}\"')
        extract_to_files = True

    # Whether to add absolute configuration assignment to the final table via getconfig

    if "c" in argv:
        show_config = True
        from getconfig import get_absolute
        argv.remove("c")
    else:
        show_config = False

    # Start extracting stuff

    jobs = []
    for name in argv[1:]:
        print(f'Reading {name}...', end="\r")
        try:
            ee = float(getoutput(f'grep \""FINAL SINGLE POINT ENERGY"\" {name} | tail -1').rstrip(" Eh").split()[-1])
            job = Job(name, free_energy=g)
            job.electronic_energy = ee

            if g:
                d = freqdir + '/' if composite else ''
                gcorr = float(getoutput(f'grep \"G-E(el)\" {d}{name} | tail -1').split()[2])
                job.gcorr = gcorr
                job.free_energy = ee + gcorr

            jobs.append(job)

        except (IndexError, ValueError):
            pass

    if not jobs:
        raise FileNotFoundError('No jobs found.')

    print()
    jobs = sorted(jobs)

    if composite:
        print(f"\nGrepped SP EE from this folder and gcorrs from {freqdir}\n")
    else:
        grepped = 'Free Energies' if g else 'Electronic Energies'
        print(f"\nGrepped {grepped}\n")

    maxlen = max([len(job.name) for job in jobs]) + 2
    min_e = min([job.get_comparison_energy() for job in jobs])

    for job in jobs:
        converged = True if getoutput(f'grep HURRAY {job.name}') != '' else False
        iterations = "conv" if converged else getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {job.name} -c')
        running = ', Running' if scratchdir in job.name else ''

        if show_config:
            try:
                coords = read_xyz(f'{job.name.split(".")[0]}.xyz').atomcoords[-1]
                job.config = get_absolute(coords, thr=60)
                
            except AssertionError:
                job.config = '?'

            config = f'- {job.config} '
        else:
            config = ''

        print(f'{job.name:<{maxlen}} - {round((job.get_comparison_energy()-min_e)*EH_TO_KCAL, 3): 10} {config}- [{iterations}{"" if converged else " iterations"}{running}]')

    ### Print table
    print()

    from prettytable import PrettyTable

    table = PrettyTable()
    table.field_names = ['#', 'Filename', 'Electronic Energy (Eh)']

    for i, job in enumerate(jobs):
        table.add_row([i+1, job.name, job.electronic_energy])

    if g:
        table.add_column('G_corr (Eh)', [job.gcorr for job in jobs])
        table.add_column('G (Eh)', [job.gcorr+job.electronic_energy for job in jobs])
        
    letter = 'G' if g else 'EE'
    table.add_column(f'Rel. {letter} (kcal/mol)', [round((job.get_comparison_energy()-min_e)*EH_TO_KCAL, 2) for job in jobs])

    if show_config:
        table.add_column('Abs. Config.', [job.config for job in jobs])

    print(table.get_string())

    ### Extract structures to file(s)

    if extract_to_files:

        def make_xyz_from_orca_inp(path):
            with open(path, 'r') as f:
                s = f.read()

            xyzblock = '\n'.join(s.split('*xyz')[1].split('\n')[1:-2])

            write_to = os.path.join(os.path.dirname(path), 'sp.xyz')
            with open(write_to, 'w') as f:
                f.write(str(len(xyzblock.split('\n'))))
                f.write('\n\n')
                f.write(xyzblock)

        import os
        import numpy as np
        from prune import cl_similarity_refining

        before = len(jobs)
        jobs = [job for job in jobs if (job.get_comparison_energy()-min_e)*EH_TO_KCAL < energy_thr]

        print(f"Keeping {len(jobs)}/{before} structures ({energy_thr} kcal/mol thr)")
        print('Reading optimized structures...')

        for job in jobs:

            trajname = job.name.split('.')[0]+"_trj.xyz"
            if os.path.isfile(trajname):
                filename = trajname
            else:
                filename = job.name.split('.')[0]+".xyz"
                if os.path.isfile(filename):
                    print(f'Trajectory file ({trajname}) not found for {job.name} - using {filename}')
                
                inp_path = os.path.join(os.path.dirname(filename), 'inp')
                if os.path.isfile(inp_path):
                    print(f'Writing sp xyzfile for {filename}...')
                    make_xyz_from_orca_inp(inp_path)

            mol = read_xyz(filename)

            job.atomnos = mol.atomnos
            job.last_coords = mol.atomcoords[-1]

        coords = np.array([job.last_coords for job in jobs])
        ref_graph = graphize(jobs[0].last_coords, jobs[0].atomnos)
        
        print('Removing similar structures...')
        coords, payload = cl_similarity_refining(coords, jobs[0].atomnos, ref_graph, payload=[np.array(jobs, dtype=object)])
        jobs = payload[0]

        # in case we discarded the most stable, recompute and save Rel. E.s
        min_e = min([job.get_comparison_energy() for job in jobs])
        for job in jobs:
            job.rel_energy = (job.get_comparison_energy()-min_e)*EH_TO_KCAL

        if os.path.basename(outname) == '':
            print(f'--> Transfering {len(jobs)} structures to same-name new files in {os.path.dirname(outname)}')

            for job in jobs:

                outname = os.path.join(os.path.dirname(outname), job.name)[:-4] + '.xyz'
                with open(outname, "w") as f:
                    title = f"Rel. {letter}. = {job.rel_energy:.3f} kcal/mol"
                    if g:
                        title += f" - gcorr(Eh) = {job.gcorr:8}"
                    write_xyz(job.last_coords, job.atomnos, f, title=title)

        else:
            with open(outname, "w") as f:
                for job in jobs:
                    title = f"{job.name} - Rel. {letter}. = {job.rel_energy:.3f} kcal/mol"
                    if g:
                        title += f" - gcorr(Eh) = {job.gcorr:8}"
                    write_xyz(job.last_coords, job.atomnos, f, title=title)

            print(f"Wrote {len(jobs)} structures to {outname}.")

if __name__ == '__main__':
    compare(sys.argv)
