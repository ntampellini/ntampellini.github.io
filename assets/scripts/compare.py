import os as operating_system
import sys
from subprocess import getoutput

import numpy as np
from firecode.pruning import get_moi_deviation_vec
from firecode.pt import pt
from firecode.rmsd import rmsd_and_max_numba
from firecode.units import EH_TO_KCAL
from firecode.utils import (graphize, read_xyz, suppress_stdout_stderr,
                            write_xyz)
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from InquirerPy.validator import PathValidator
from rich.traceback import install

install(show_locals=True, locals_max_length=None, locals_max_string=None, width=120)

class tcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[38;2;65;165;165m'
    WARNING = '\033[93m'
    FAIL = '\033[38;2;232;111;136m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

scratchdir = '/vast/palmer/scratch/miller/nt383'

class Options:
    def __init__(self):
        return

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
    
    def get_coords(self):
        if not hasattr(self, "coords"):
            self.read_coords()

        return self.last_coords
        
    def read_coords(self):
        if self.name.endswith(".out"):
            self.last_coords = read_xyz(f'{self.name[:-4]}.xyz').atomcoords[-1]

        elif self.name.endswith(".xyz"):
            raise NotImplementedError

def compare(argv):
    '''
    '''

    # clean input files
    for filename in argv[1:]:
        if not (filename.endswith(".out") or filename.endswith(".xyz")): 
            argv.remove(filename)

    if len(argv) == 1:
        print("\n  Compare completed ORCA jobs electronic/free energy of groups of jobs and prints " +
            "collected energy values. Extra extraction options. Syntax:\n\n" +
            "  python compare.py *.out")
        quit()

    else:
        print(f"--> Specified {len(argv)} outfiles.")

    ### Multiple option selector
    options = Options()
    
    # check if we have G(corr) values to grep before the user has a chance to choose free energy
    avail_gcorrs = []
    if "GIBBS" in getoutput(f'grep GIBBS {argv[1]}'):
        avail_gcorrs.append(Choice(value=".", name=f'./    This folder   - read free energy from this folder.'))

    if "GIBBS" in getoutput(f'grep GIBBS ../{argv[1]}'):
        avail_gcorrs.append(Choice(value="..", name=f'../   Parent folder - read free energy from the parent folder.'))

    avail_gcorrs.append(Choice(value=None,    name=f'?     Other folder  - choose another folder to read G(corr) values.'))

    choices = [
        Choice(value="g",name='G          - Use free energy as a comparison metric.', enabled=len(avail_gcorrs)>1),
        Choice(value="x",          name='extract    - extract structures to one or more new files.', enabled=False),
        Choice(value="stereochem", name='stereochem - calculate/show stereochemistry around a set of atoms.', enabled=False),
        Choice(value="novel",      name='novel      - only consider structures that are dissimilar to some others.', enabled=False),
        Choice(value="scratch",    name='scratch    - scrape data from currently running jobs that were started in this folder.', enabled=False),
    ]

    options_to_set = inquirer.checkbox(
            message="Select options (spacebar to toggle, enter to confirm):",
            choices=choices,
            cycle=False,
            disabled_symbol='⬡',
            enabled_symbol='⬢',
        ).execute()
    
    for choice in choices:
        setattr(options, choice.value, choice.value in options_to_set)

    ### Read structures to compare to the specified ones,
    ### to filter the ones already present
    if options.novel:

        comparison_structs = []

        while True:
            filename = inquirer.filepath(
                message="Select a structure file for comparison:",
                default="./" if operating_system.name == "posix" else "C:\\",
                validate=PathValidator(is_file=True, message="Input is not a file"),
                only_files=True,
            ).execute()

            comparison_structs.extend(read_xyz(filename).atomcoords)

            if not inquirer.confirm(
                message="Specify more files?",
                default=False).execute():
                break

    ### If checking scratch, append the running job full names to argv
    running_names = []
    if options.scratch:
       
        cwd = operating_system.getcwd()

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
                    running_names.append(fullname)

                    with suppress_stdout_stderr():
                        # copy traj from scratch
                        operating_system.system(f'cp {scratchdir}/{jobdir}/{jobname}_trj.xyz .')
                        operating_system.system(f'cp {scratchdir}/{jobdir}/{jobname}_Compound_*_trj.xyz .')
                        operating_system.system(f'cp {scratchdir}/{jobdir}/{jobname}.property.txt .')
                        # operating_system.system(f'cp {scratchdir}/{jobdir}/{jobname}.out .')

            except StopIteration:
                continue

    ### If interested in a composite method, see if the
    ### user specified where to extract gcorrs

    # if "freq" in [kw.split("=")[0] for kw in argv]:
    #     freqdir = next((kw.split("=")[-1] for kw in argv if "freq" in kw))
    #     argv.remove(f"freq={freqdir}")
    # else:
    #     freqdir = '..'

    ### If we are interested in free energy, set appropriate lookup folder
    if options.g:

        freqdir = inquirer.select(
            message="Which folder would you like to extract G(corr) values from?",
            choices=avail_gcorrs,
            default=avail_gcorrs[0].value,
        ).execute()

        if freqdir is None:
            freqdir = inquirer.filepath(
                message="Pick a folder to extract G(corr) values from:",
                only_directories=True,
                validate=PathValidator(is_dir=True),
        ).execute()

        # folderstring = 'parent folder(s)' if freqdir == '..' else freqdir
        # print(f'--> Extracting G(corr) values from {folderstring}')

    ### Set the energy threshold for file extraction

    # if "thr" in [kw.split("=")[0] for kw in argv]:
    #     energy_thr = next((kw.split("=")[-1] for kw in argv if "thr=" in kw))
    #     print(f'-> set thr to {energy_thr} kcal/mol')
    #     try:
    #         argv.remove(f"thr={energy_thr}")
    #     except ValueError:
    #         raise ValueError(f'Trying to delete non-existent\"thr={energy_thr}\"')

    if options.x:
        energy_thr = inquirer.text(
            message="Extraction threshold from the most stable? (kcal/mol)",
            default="5",
            validate=lambda x: x.replace(".", "").isdigit(),
            filter=lambda x: float(x),
        ).execute()

    # Specify where to extract files - if a filename is provided, structures
    # will be extracted to a single file, otherwise the same filenames of
    # initial structures will be used

    # if "x" in [kw.split("=")[0] for kw in argv]:
    #     outname = next((kw.split("=")[-1] for kw in argv if "=" in kw))
    #     try:
    #         argv.remove(f"x={outname}")
    #     except ValueError:
    #         raise ValueError(f'Trying to delete non-existent\"x={outname}\"')
    #     extract_to_files = True

        outfolder = inquirer.filepath(
            message="Where do you want to extract structures?",
            default=operating_system.getcwd(),
            validate=PathValidator(is_dir=True, message="Please specify a directory."),
            only_directories=True,
        ).execute()

        outfile_name = inquirer.text(
            message="Filename to save structures to? (leave blank for same-name multiple files)",
            default="",
        ).execute()

        outname = operating_system.path.join(outfolder, outfile_name)

    # Whether to add absolute configuration assignment to the final table via getconfig

    if options.stereochem:
        from getconfig import get_absolute

    # check if we are reading compound jobs,
    # and if so ask what energy to extract
    previous_to_last_energy = False
    try:
        if operating_system.path.isfile(f'{argv[1][:-4]}.inp'):
            if getoutput(f'grep \".cmp\" {argv[1][:-4]}.inp'):
                previous_to_last_energy = inquirer.select(
                    message="What energy would you like to extract?",
                    choices=(
                        Choice(value=True, name='Previous to last energy (last at lower level)'),
                        Choice(value=False, name='Last energy (higher level)'),
                    ),
                    default=False,
                ).execute()

                freqdir = "."
                
    except IndexError:
        pass

    # Start extracting stuff
    jobs = []
    for name in argv[1:]:
        print(f'Reading {name}...', end="\r")

        # parse energies from .xyz comment line
        if name.endswith(".xyz"):

            # get lines right after the number of atom, which should contain the energy
            energies = getoutput(f'grep -A1 "^[0-9]\+$" {name} | grep -v "^[0-9]\+$" | grep -v "^--$"').split("\n")

            relative = any(['Rel' in line for line in energies])

            energies = [float(next(
                            part for part in line.split() if (any(c == '.' for c in part) and
                                                              set(part).issubset('0123456789.-'))))
                        for line in energies]
            
            # with .xyz files we read the structures here so that we can deal with
            # multiple structures in a single file
            mol = read_xyz(name)
            assert len(mol.atomcoords) == len(energies)
            
            # if they are Relative, chances are they are in kcal/mol
            if relative:
                energies = [e/EH_TO_KCAL for e in energies]

            for i, energy in enumerate(energies):
                job = Job(name[:-4]+f"_conf{i}", free_energy=options.g)
                job.electronic_energy = energy
                job.last_coords = mol.atomcoords[i]
                job.atomnos = mol.atomnos
                jobs.append(job)

        # parse ORCA output files 
        elif name.endswith(".out"):
            try:

                if not previous_to_last_energy:
                    ee = float(getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {name} | tail -1').rstrip(" Eh").split()[-1])

                else:
                    ee = float(getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {name} | tail -2 | head -1').rstrip(" Eh").split()[-1])

                job = Job(name, free_energy=options.g)
                job.electronic_energy = ee

                if options.g:
                    # d = freqdir + '/' if composite else ''
                    current = operating_system.getcwd()
                    file_folder = operating_system.path.dirname(name)
                    if file_folder:
                        operating_system.chdir(file_folder)
                    operating_system.chdir(freqdir)
                    target_folder = operating_system.getcwd()
                    target_filename = operating_system.path.basename(name)
                    target_filename_abs = operating_system.path.join(target_folder, target_filename)
                    operating_system.chdir(current)

                    gcorr = float(getoutput(f'grep \"G-E(el)\" {target_filename_abs} | tail -1').split()[2])
                    job.gcorr = gcorr
                    job.free_energy = ee + gcorr

                jobs.append(job)

            except (IndexError, ValueError) as e:
                print(e)
                pass

    if not jobs:
        raise FileNotFoundError('No jobs found.')

    print()
    jobs = sorted(jobs)

    if options.g:
        print(f"\nGrepped SP EE from this folder and G(corr) values from {freqdir}\n")

    maxlen = max([len(job.name) for job in jobs]) + 2
    min_e = min([job.get_comparison_energy() for job in jobs])

    for job in jobs:
        converged = getoutput(f'grep HURRAY {job.name}') != ''
        iterations = "conv" if converged else getoutput(f'grep \"FINAL SINGLE POINT ENERGY\" {job.name} -c')
        running = ', Running' if scratchdir in job.name else ''

        # save job.completed attribute
        look_for = "ORCA TERMINATED NORMALLY"
        job.completed = getoutput(f'grep \'{look_for}\' {job.name}') != ''

        if options.stereochem:
            try:
                coords = job.get_coords()
                job.config = get_absolute(coords, thr=60)
                
            except AssertionError:
                job.config = '?'

            config = f'- {job.config} '
        else:
            config = ''

        if options.novel:

            last_coords = job.get_coords()
            masses = np.array([pt[a].mass for a in job.atomnos])
            for reference in comparison_structs:
                moi_dev_vec = get_moi_deviation_vec(last_coords, reference, masses)
                if np.all(moi_dev_vec < 0.01):
                    job.novelty = False
                    break

                rmsd, maxdev = rmsd_and_max_numba(last_coords, reference, center=True)
                if rmsd < 0.25 and maxdev < 0.5:
                    job.novelty = False
                    break

                job.novelty = True

        print(f'{job.name:<{maxlen}} - {round((job.get_comparison_energy()-min_e)*EH_TO_KCAL, 3): 10} {config}- [{iterations}{"" if converged else " iterations"}{running}]')

    ### Print table
    print()

    from prettytable import PrettyTable

    table = PrettyTable()
    table.field_names = ['#', 'Filename', 'Electronic Energy (Eh)']

    for i, job in enumerate(jobs):
        table.add_row([i+1, job.name, job.electronic_energy])

    if options.g:
        table.add_column('G_corr (Eh)', [job.gcorr for job in jobs])
        table.add_column('G (Eh)', [job.gcorr+job.electronic_energy for job in jobs])
        
    letter = 'G' if options.g else 'EE'
    table.add_column(f'Rel. {letter} (kcal/mol)', [round((job.get_comparison_energy()-min_e)*EH_TO_KCAL, 2) for job in jobs])
    table.add_column('Status', [tcolors.OKGREEN + tcolors.BOLD + "COMPLETED" + tcolors.ENDC if job.completed
                           else (tcolors.WARNING + "RUNNING" + tcolors.ENDC if job.name in running_names
                            else tcolors.FAIL + "INCOMPLETE" + tcolors.ENDC) for job in jobs])

    if options.stereochem:
        table.add_column('Abs. Config.', [job.config for job in jobs])

    if options.novel:
        table.add_column('Novelty', ["NEW" if job.novelty else "   " for job in jobs])

    print(table.get_string())

    if options.stereochem:
        major_config = min(jobs, key=lambda job: job.get_comparison_energy()).config
        minor_config = 'R' if major_config == 'S' else 'S'
        configs = [job.config for job in jobs]
        count = configs.count(major_config)

        major_energies = [job.get_comparison_energy() for job in jobs if job.config == major_config]
        minor_energies = [job.get_comparison_energy() for job in jobs if job.config == minor_config]

        major_lowest_e = min(major_energies)
        minor_lowest_e = min(minor_energies)
        delta_kcal = (minor_lowest_e-major_lowest_e)*EH_TO_KCAL

        R = 0.001985877534
        T = (273.15+25)

        k = np.exp(-delta_kcal/(R*T))
        major_amt = round(100 / (k+1), 1)
        minor_amt = round(100.0 - major_amt, 1)

        # now again but cumulative on all confs
        rel_e_major_config = [(e-major_lowest_e)*EH_TO_KCAL for e in major_energies]
        pop_major = sum([np.exp(-e/(R*T)) for e in rel_e_major_config])

        rel_e_minor_config = [(e-major_lowest_e)*EH_TO_KCAL for e in minor_energies]
        pop_minor = sum([np.exp(-e/(R*T)) for e in rel_e_minor_config])

        major_amt_cumulative = round(pop_major/(pop_major+pop_minor)*100, 1)
        minor_amt_cumulative = round(100.0 - major_amt_cumulative, 1)

        print(f'Major config: {major_config} ({count}/{len(jobs)})')
        print(f'  pred. {major_amt:>4}:{minor_amt:>4} er (from the two best states)')
        print(f'  pred. {major_amt_cumulative:>4}:{minor_amt_cumulative:>4} er (from all the states)')

    ### Extract structures to file(s)

    if options.x:

        def make_xyz_from_orca_inp(path):
            with open(path, 'r') as f:
                s = f.read()

            xyzblock = '\n'.join(s.split('*xyz')[1].split('\n')[1:-2])

            write_to = operating_system.path.join(operating_system.path.dirname(path), 'sp.xyz')
            with open(write_to, 'w') as f:
                f.write(str(len(xyzblock.split('\n'))))
                f.write('\n\n')
                f.write(xyzblock)

        from prune import cl_similarity_refining

        before = len(jobs)
        jobs = [job for job in jobs if (job.get_comparison_energy()-min_e)*EH_TO_KCAL < energy_thr]

        print(f"Keeping {len(jobs)}/{before} structures ({energy_thr} kcal/mol thr)")
        print('Reading optimized structures...')

        for job in jobs:
            if not hasattr(job, "last_coords"):

                trajname = job.name[:-4]+"_trj.xyz"
                if operating_system.path.isfile(trajname):
                    filename = trajname

                elif operating_system.path.isfile(trajname := job.name[:-4]+"_Compound_1_trj.xyz"):
                    filename = trajname

                else:
                    filename = job.name[:-4]+".xyz"
                    if operating_system.path.isfile(filename):
                        print(f'Trajectory file ({trajname}) not found for {job.name} - using {filename}')
                    
                    inp_path = operating_system.path.join(operating_system.path.dirname(filename), 'inp')
                    if operating_system.path.isfile(inp_path):
                        print(f'Writing sp xyzfile for {filename}...')
                        make_xyz_from_orca_inp(inp_path)

                mol = read_xyz(filename)

                job.atomnos = mol.atomnos
                job.last_coords = mol.atomcoords[-1]
        
        if options.novel:
            before = len(jobs)
            jobs = [job for job in jobs if job.novelty]
            print(f'Considering only novel structures: retaining {len(jobs)}/{before}')

        coords = np.array([job.last_coords for job in jobs])
        ref_graph = graphize(jobs[0].last_coords, jobs[0].atomnos)

        print('Removing similar structures...')
        coords, payload = cl_similarity_refining(coords, jobs[0].atomnos, ref_graph, payload=[np.array(jobs, dtype=object)])
        jobs = payload[0]

        # in case we discarded the most stable, recompute and save Rel. E.s
        min_e = min([job.get_comparison_energy() for job in jobs])
        for job in jobs:
            job.rel_energy = (job.get_comparison_energy()-min_e)*EH_TO_KCAL

        # print to same-name new files if user asked to
        if operating_system.path.basename(outname) == '':
            print(f'--> Transfering {len(jobs)} structures to same-name new files in {operating_system.path.dirname(outname)}')

            for job in jobs:

                outname = operating_system.path.join(operating_system.path.dirname(outname), job.name).rstrip('.xyz').rstrip('.out') + '.xyz'
                with open(outname, "w") as f:
                    title = f"{letter} = {job.get_comparison_energy()} Eh - Rel. {letter}. = {job.rel_energy:.3f} kcal/mol"
                    if options.g:
                        title += f" - gcorr(Eh) = {job.gcorr:8}"
                    write_xyz(job.last_coords, job.atomnos, f, title=title)

        # printing to a single new filename the user provided
        else:
            with open(outname, "w") as f:
                for job in jobs:
                    title = f"{job.name} - {letter} = {job.get_comparison_energy()} Eh - Rel. {letter}. = {job.rel_energy:.3f} kcal/mol"
                    if options.g:
                        title += f" - G(corr)(Eh) = {job.gcorr:8}"
                    write_xyz(job.last_coords, job.atomnos, f, title=title)

            print(f"Wrote {len(jobs)} structures to {outname}.")

if __name__ == '__main__':
    compare(sys.argv)
