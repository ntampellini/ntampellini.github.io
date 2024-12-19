####################################################
hourly_rate = 0.004 # $/(h*core)

options = {
    "solvent" : "chloroform",
    "solvent_model" : 'CPCM',
    "level" : "R2SCAN-3c",
    "basis_set" : "",
    "opt": "",
    "freq" : False,
    "temp" : 273.15+25,
    "procs" : 16,
    "mem" : 6, # Memory per core, in GB
    "charge" : 0,
    "maxstep" : 0.05, # in Bohr atomic units (1au = 0.529177 A)
    "popt" : False,
    "compound" : False,
    "compound_job_scriptname" : "optf+sp.cmp",
    "hh" : False, # hybrid hessian

    "additional_kw" : "Defgrid3",
    "extra_block" : "",
    }

# pick name from ORCA 6.0.0 manual (p. 1006) and
# epsilon value from https://people.chem.umass.edu/xray/solvent.html
epsilon_dict = {
    "acetonitrile" : 37.5,
    "benzene" : 2.3,
    "chloroform" : 4.81,
    "dichloromethane" : 9.04,
    "diethyl ether" : 4.33,
    "diethylether" : 4.33,
    "dimethylformamide" : 36.71,
    "dimethylacetamide" : 37.78,
    "thf" : 7.58,
    "tetrahydrofuran" : 7.58,
    "dmf" : 36.71,
    "dmso" : 47.2,
    "ethyl acetate" : 6.02,
    "ethyl ethanoate" : 6.02,
    "phcf3" : 9.18,
    "toluene" : 2.38,
}

# round temperature so it looks prettier
options["temp"] = round(options["temp"], 2)

####################################################

import argparse
import os
import re
import shutil
import sys
import numpy as np
from subprocess import getoutput, run

from InquirerPy import inquirer
from InquirerPy.base.control import Choice

from utils import norm_of, pt, read_xyz, suppress_stdout_stderr, point_angle, dihedral

def inquirer_set_options(args):
    '''
    InquirerPy command-line option setter
    
    '''
    print()

    runtype = inquirer.select(
        message="Which kind of input file would you like to generate?",
        choices=(
            Choice(value='fastsp',   name='fastsp   - Low-level DFT single-point energy calculation.'),
            Choice(value='sp',       name='sp       - High-level DFT single-point energy calculation.'),
            Choice(value='optf',     name='optf     - Geom. optimization + frequency calculation.'),
            Choice(value='popt',     name='popt     - Partial optimization (specify constraints).'),
            Choice(value='ts',       name='ts       - Saddle optimization + frequency calculation.'),
            Choice(value='scan',     name='scan     - Perform a distance/angle/dihedral scan.'),
            Choice(value='nmr',      name='nmr      - Single-point NMR tensors calculation.'),
            Choice(value='compound', name='compound - Choose a compound method routine.'),
            Choice(value='irc',      name='irc      - Intrinsic reaction coordinate calculation.'),
            Choice(value='tddft',    name='tddft    - TD-DFT calculation.'),
        ),
        default='sp',
    ).execute()

    # modify the option on the args namespace
    setattr(args, runtype, True)

    # set or confirm solvent
    options["solvent"] = inquirer.fuzzy(
        message=f"Current solvent is {options['solvent']}. Press enter to confirm or type another solvent:",
        choices=epsilon_dict.keys(),
        default=options["solvent"],
        # validate=lambda result: result in epsilon_dict,
        # invalid_message="Solvent not recognized."
    ).execute()

def get_comp_script_inp(rootname):
    return f'''* xyzfile {options["charge"]} 1 {rootname}.xyz
    
%pal
 nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

%geom
  MaxStep {options["maxstep"]}
end

%compound "optf+sp.cmp"
  with
    solvent = "{options["solvent"]}";
end

'''

def get_cpcm_block(solvent, solvent_model):
    s = '%cpcm\n'
    if solvent_model == 'CPCM':
        s += f"  epsilon {epsilon_dict[options['solvent']]}\n"
    elif solvent_model == 'SMD':
        s += f"  smd True\n  SMDsolvent \"{options['solvent']}\"\n"
    else:
        raise Exception(f'Solvent model \"{solvent_model}\" not recognized (use CPCM or SMD)')

    s += "end"

    return s

def get_inp(rootname, opt=True):
    maxstep_line = f'MaxStep {options["maxstep"]}'
    hess_line = "Calc_Hess true" if args.ts else ""

    if options["hh"]:
        hess_line += f'\n  Hybrid_Hess {"{"+hh_ids+"}"} end'

    geom_block = f'%geom\n\n  {maxstep_line}\n\n  {hess_line}\n\nend' if opt else ""

    inp = f'''! {options["level"]} {options["basis_set"]} {"CPCM" if options["solvent_model"] == "CPCM" else ""} {options["opt"]}
! {options["additional_kw"]}

%pal
  nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

{geom_block}

{get_cpcm_block(options["solvent"], options["solvent_model"])}

{f"%freq temp {options['temp']} end" if options["freq"] else ""}

{options["extra_block"]}

* xyzfile {options["charge"]} 1 {rootname}.xyz

'''

    return inp

def get_scan_block(filename):

    # get scan type
    c_type = inquirer.select(
            message="Which type of scan would you like to perform?",
            choices=(
                Choice(value='B', name='Bond (2 indices, target distance)'),
                Choice(value='A', name='Angle (3 indices, target angle)'),
                Choice(value='D', name='Dihedral (4 indices, target dihedral)'),
            ),
            default='B',
        ).execute()

    required_indices = 2 if c_type == "B" else (3 if c_type == "A" else 4)

    # get scan indices and target
    while True:
        target_modifier_string = ' (or put \'ts\' as target)' if c_type == 'B' else (' (or leave target blank for a 360° rotation)' if c_type == 'D' else '')
        c_string = inquirer.text(
                message=f"Specify {required_indices} indices and a target{target_modifier_string}:",
            ).execute()
        
        if c_type == 'D' and  len(c_string.split()) == 4:
            c_string += ' full'

        if len(c_string.split()) == required_indices+1:
            indices = [int(n) for n in c_string.split()[:-1]]
            break        
    
    # read coordinates
    mol = read_xyz(filename)
    coords = mol.atomcoords[-1]

    # get start and end values, dertermine number of steps
    if c_type == "B":
        i1, i2 = indices
        e1, e2 = pt[mol.atomnos[i1]], pt[mol.atomnos[i2]]
        start_value = round(norm_of(coords[i1]-coords[i2]), 2)

        if c_string.split()[-1] == "ts":
            ts_d_estimate = get_ts_d_estimate(filename, indices)
            sign = np.sign(ts_d_estimate - start_value)
            target_value = round(ts_d_estimate + 0.2 * sign, 2)

            elem_id_string = f'{e1}{i1}-{e2}{i2}'

        else:
            target_value = round(float(c_string.split()[-1]))

    if c_type == "A":
        i1, i2, i3 = indices
        e1, e2, e3 = pt[mol.atomnos[i1]], pt[mol.atomnos[i2]], pt[mol.atomnos[i3]]
        start_value = round(point_angle(coords[i1], coords[i2], coords[i3]), 2)
        target_value = round(float(c_string.split()[-1]))

        elem_id_string = f'{e1}{i1}-{e2}{i2}-{e3}{i3}'

    if c_type == "D":
        i1, i2, i3, i4 = indices
        e1, e2, e3, e4 = pt[mol.atomnos[i1]], pt[mol.atomnos[i2]], pt[mol.atomnos[i3]], pt[mol.atomnos[i4]]
        start_value = round(dihedral(coords[np.array(indices)]), 2)

        if c_string.split()[-1] == 'full':
            target_value = round(start_value + 360)
        else:
            target_value = round(float(c_string.split()[-1]))

        elem_id_string = f'{e1}{i1}-{e2}{i2}-{e3}{i3}-{e4}{i4}'

    stepsize = 0.05 if c_type == 'B' else (5 if c_type == 'A' else 10)

    n_steps = round((target_value-start_value)/stepsize)
    new_target = round(start_value + n_steps * stepsize, 2)
    n_steps = abs(n_steps)+1

    output_string = '%geom Scan\n'
    output_string += f'        {c_type} {" ".join([str(i) for i in indices])} = {start_value}, {new_target}, {n_steps}\n'
    output_string += '        end\n      end\n'

    uom = ' Å' if c_type == 'B' else '°'

    print(f'--> Scanning {elem_id_string} from {start_value}{uom} to {target_value}{uom} in {n_steps} steps.')

    return output_string

def get_procs(rootname):
    '''
    
    '''
    natoms = float(getoutput(f'head {rootname}.xyz').split('\n')[0])
    
    if natoms < 100:
        procs = 4

    else:
        procs = 8

    if options["opt"]:
        procs *= 2

    if options["freq"]:
        procs *= 2

    return min(procs, 32)

def get_daily_cost(procs, mem_per_core):
    if not args.priority:
        return 0
    else:
        return max(procs, procs*mem_per_core/15) * 24 * hourly_rate

def get_ts_d_estimate(filename, indices, factor=1.35, verbose=True):
    '''
    Returns an estimate for the distance between two
    specific atoms in a transition state, by multipling
    the sum of covalent radii for a constant.
    
    '''
    mol = read_xyz(filename)
    i1, i2 = indices
    a1, a2 = pt[mol.atomnos[i1]], pt[mol.atomnos[i2]]
    cr1 = a1.covalent_radius
    cr2 = a2.covalent_radius

    est_d = round(factor * (cr1 + cr2), 2)

    if verbose:
        print(f'--> Estimated TS d({a1}-{a2}) = {est_d} Å')
        
    return est_d

############################################################# START OF LOGIC

# usage = (
#     "\n  Makes one or more ORCA inputs following the desired specifications. Syntax:\n\n" +
#     "  python mkorca.py conf*.xyz [option]\n\n" + 
#     "  conf*.xyz: base name of input geometry file(s)\n" +
#     "  option:\n" +
#     "    sp:     single-point energy calculation.\n" +
#     "    fastsp: fast single-point energy calculation.\n" +
#     "    optf:   optimization + frequency calculation.\n" +
#     "    popt:   partial optimization (specify distance/angle/dihedral constraints).\n" +
#     "    ts:     eigenvector-following saddle point optimization.\n" +
#     "    nmr:    single-point energy calculation with NMR shieldings.\n" + 
#     "  Each of these options might have different default levels of theory. Manually check/modify this script at your convenience.\n"
# )

parser = argparse.ArgumentParser() 
parser.add_argument("inputfiles", help="Input filenames, in .xyz format.", action='store', nargs='*', default=None)
parser.add_argument("-solvent", help="Set solvent to the specified one.", action="store", required=False)
parser.add_argument("-sp", help="Set input type to SP.", action="store_true", required=False)
parser.add_argument("-fastsp", help="Set input type to fast SP.", action="store_true", required=False)
parser.add_argument("-ts", help="Set input type to TS.", action="store_true", required=False)
parser.add_argument("-hh", help="Specify the use of a hybrid hessian.", action="store_true", required=False)
parser.add_argument("-optf", help="Set input type to optimization + frequency calculation.", action="store_true", required=False)
parser.add_argument("-popt",help="Set input type to partial optimization.", action="store_true", required=False)
parser.add_argument("-scan", help="Perform a distance/angle/dihedral scan.", action="store_true", required=False)
parser.add_argument("-nmr", help="Set input type to NMR (single-point).", action="store_true", required=False)
parser.add_argument("-compound", help="Set input type to a compound method.", action="store_true", required=False)
parser.add_argument("-irc", help="Set input type to IRC.", action="store_true", required=False)
parser.add_argument("-priority", help="Run jobs with priority (and estimate cost)", action="store_true", required=False)
args = parser.parse_args()

if not args.inputfiles:
    raise Exception('No input structures selected. Please specify at least one.')

if not any((args.sp, args.fastsp, args.ts, args.optf, args.popt, args.scan, args.nmr, args.compound, args.irc)):
    inquirer_set_options(args)

if args.solvent:
    print(f"--> Setting solvent to {options['solvent']}")

if args.sp:
    options["opt"] = ""
    options["freq"] = False
    options["level"] = 'wB97M-V'
    options["basis_set"] = 'def2-TZVPP'
    options["solvent_model"] = 'SMD'

if args.fastsp:
    options["opt"] = ""
    options["freq"] = False
    options["level"] = 'R2SCAN-3c'
    options["basis_set"] = ''

if args.ts:
    options["opt"] = "OptTS"
    options["mem"] = 8
    options["freq"] = True
    options["additional_kw"] += " TightOpt LARGEPRINT"

if args.hh:
   
    hh_ids = input("Hybrid Hessian: provide indices involved in TS: ")

    if hh_ids == "":
        popt_inp_name = args.inputfiles[0].split('.')[0] + '.inp'
        if popt_inp_name in os.listdir(os.path.dirname(os.getcwd())):
            popt_ids = getoutput(f'grep {{*}} ../{popt_inp_name}')
            popt_ids = ' '.join(re.findall(r'(?<!\.)\b\d+\b(?!\.)', popt_ids))
            popt_ids = ' '.join({i for i in popt_ids.split()})
            ans = input(f"Would you like to use indices from the previous popt? [{popt_ids}] [y]/n")
            if ans in ('y', 'Y', ''):
                hh_ids = popt_ids

            else:
                print('Please specify hybrid hessian-excluded indices.')
                sys.exit()
        else:
            print('Please specify hybrid hessian-excluded indices.')
            sys.exit()

if args.optf:
    options["opt"] = "Opt"
    options["freq"] = True
    options["additional_kw"] += " TightOpt LARGEPRINT"

if args.popt:
    options["freq"] = False
    options["opt"] = "Opt"
    options["popt"] = True
    # options["additional_kw"] += " TightOpt LARGEPRINT"

    constraint_strings = []
    while True:

        c_type = inquirer.select(
            message="Select which type of constraint you would like to specify:",
            choices=(
                Choice(value=None, name='(End of constraints)'),
                Choice(value='B', name='Bond (2 indices, optional distance)'),
                Choice(value='A', name='Angle (3 indices, optional angle)'),
                Choice(value='D', name='Dihedral (4 indices, optional dihedral)'),
            ),
            default='B',
        ).execute()

        if c_type is None:
            break

        autod_ts_line = ", or \"ts\" to estimate one" if c_type == "B" else ""

        c_string = inquirer.text(
            message=f"Provide indices and optional distance{autod_ts_line}:"
            ).execute()

        if c_type == "B" and c_string.split()[-1] == "ts":
            indices = [int(i) for i in c_string.split()[:-1]]
            ts_d_estimate = str(get_ts_d_estimate(args.inputfiles[0], indices))
            c_string = ' '.join(c_string.split()[:-1] + [ts_d_estimate])

        constraint_strings.append((c_type, c_string))

    options["extra_block"] = "%geom Constraints\n"
    for c_type, c_string in constraint_strings:
        options["extra_block"] += "  {{ {0} {1} C }}\n".format(c_type, c_string)
        
    options["extra_block"] += "  end\nend"

if args.scan:
    options["opt"] = "Opt"

    # make scan %geom block and add it to the extra block
    options["extra_block"] += get_scan_block(args.inputfiles[0])

if args.nmr:
    options["freq"] = False
    options["opt"] = ""
    options["level"] = 'PBE0'
    options["basis_set"] = '6-311+G(2d,p)'
    options["additional_kw"] += " NMR"

if args.compound:
    try:
        script_folder = os.path.dirname(os.path.realpath(__file__))
        source = os.path.join(script_folder, options["compound_job_scriptname"])
        shutil.copyfile(source, options["compound_job_scriptname"])
        print(f'--> Copied {options["compound_job_scriptname"]} to input file folder.')
    except Exception:
        raise Exception(f'Something went wrong when copying {options["compound_job_scriptname"]} to the destination folder.')

if args.irc:
    options["freq"] = True
    options["additional_kw"] += " IRC"
    options["opt"] = ""

if options["freq"] and "Freq" not in options["additional_kw"]:
    options["additional_kw"] += " Freq"

print()
for option, value in options.items():
   if value != "":
       print(f"--> {option} = {value}")
print()

allchars = ''.join(args.inputfiles)
auto_charges = False
if '+' in allchars or '-' in allchars:
    auto_charges = inquirer.confirm(
        message="Found charge signs in filenames (+/-). Auto assign charges in input files?",
        default=True,
    ).execute()

procs_list = [get_procs(filename.split('.')[0]) for filename in args.inputfiles]
cum_cost = 0

print(f'Filename                       Cores  Mem(GB)  Max cost (24h)')
print('-------------------------------------------------------------')

for filename, procs in zip(args.inputfiles, procs_list):
    rootname = filename.split('.')[0]
    options["procs"] = procs
    cost = get_daily_cost(procs, mem_per_core=options["mem"])
    cum_cost += cost

    print(f'{filename:30s} {procs:2}     {options["mem"]*procs:d}       {cost:.2f} $')

    if auto_charges:
        if '+' in rootname:
            options['charge'] = 1
            print(f'--> {filename} : assigned charge +1 based on input name')
        elif '-' in rootname:
            options['charge'] = -1
            print(f'--> {filename} : assigned charge -1 based on input name')

    if args.compound:
        s = get_comp_script_inp(rootname)

    else:
        s = get_inp(rootname, opt=options["opt"] != "")

    with open(f'{rootname}.inp', 'w') as f:
        f.write(s)

    # Convert all text files to Linux format
    try:
        with suppress_stdout_stderr():
            run(f'dos2unix {rootname}.inp'.split())
            run(f'dos2unix {rootname}.xyz'.split())
    except FileNotFoundError:
        pass

print('-------------------------------------------------------------')
print(f'Maximum estimated cost (24 h runtime): {cum_cost:.2f} $\n')

run_jobs = inquirer.confirm(
        message="Run jobs for the newly generated input files?",
        default=True,
    ).execute()

if run_jobs:
    inpnames = [name.split(".")[0]+".inp" for name in args.inputfiles]

    from orcasub_batch import main

    # adding dummy first element to simulate sys.argv
    main([None]+inpnames)