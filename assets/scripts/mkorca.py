####################################################
hourly_rate = 0.004 # $/(h*core)

options = {
    "solvent" : "chloroform",
    "solvent_model" : 'CPCM', # CPCM, SMD, ALPB
    "level" : "R2SCAN-3c",
    "basis_set" : "",
    "opt": "",
    "freq" : False,
    "temp" : 273.15+25,
    "procs" : 16,
    "mem" : 6, # Memory per core, in GB
    "charge" : 0,
    "mult" : 1,
    "maxstep" : 0.3, # in au, i.e. Bohr atomic units (1au = 0.529177 A) - default is 0.3 au
    "popt" : False,
    "compound" : False,
    "compound_job_scriptname" : "optf+sp.cmp",
    "comp_job_extra_variables" : "",
    "hh" : False, # hybrid hessian

    "additional_kw" : "",
    "constr_block" : "",
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
global_constraints = []

####################################################

import argparse
import os
import re
import shutil
import sys
from subprocess import getoutput, run

import numpy as np
from InquirerPy import inquirer
from InquirerPy.base.control import Choice
from rich.traceback import install

from utils import (dihedral, norm_of, point_angle, pt, read_xyz,
                   suppress_stdout_stderr)

install(show_locals=True, locals_max_length=None, locals_max_string=None, width=120)

class Constraint:
    '''
    Constraint class with indices, type and value attributes.
    
    '''
    def __init__(self, indices, value=None):
        
        self.indices = indices

        self.type = {
            2 : 'B',
            3 : 'A',
            4 : 'D',
        }[len(indices)]

        self.value = value

def inquirer_set_options(args):
    '''
    InquirerPy command-line option setter
    
    '''
    print()

    runtype = inquirer.select(
        message="Which kind of input file would you like to generate?",
        choices=(
            Choice(value='sp',       name='sp       - High-level DFT single-point energy calculation.'),
            Choice(value='optf',     name='opt(f)   - Geom. optimization (+ frequency calculation).'),
            Choice(value='popt',     name='popt     - Partial optimization (specify constraints).'),
            Choice(value='ts',       name='ts       - Saddle optimization + frequency calculation.'),
            Choice(value='goat',     name='goat     - Conformational search via GOAT.'),
            Choice(value='compound', name='compound - Choose a compound method routine.'),
            Choice(value='scan',     name='scan     - Perform a distance/angle/dihedral scan.'),
            Choice(value='fastsp',   name='fastsp   - Low-level DFT single-point energy calculation.'),
            Choice(value='nmr',      name='nmr      - Single-point NMR tensors calculation.'),
            Choice(value='irc',      name='irc      - Intrinsic reaction coordinate calculation.'),
            Choice(value='tddft',    name='tddft    - TD-DFT calculation.'),
            Choice(value='freqtemp', name='freqtemp - Recalculate vibrational corrections at a new temperature.'),
        ),
        default='sp',
    ).execute()

    # modify the option on the args namespace
    setattr(args, runtype, True)

    if not any((args.compound, args.freqtemp)):
        # set solvent model
        options["solvent_model"] = inquirer.select(
            message=f"What solvation model would you like to use?",
            choices=(
                Choice(value=None, name='vacuum'),
                Choice(value='CPCM', name='CPCM'),
                Choice(value='SMD', name='SMD'),
            ),
            default='SMD' if args.sp else 'CPCM',
            # validate=lambda result: result in epsilon_dict,
            # invalid_message="Solvent not recognized."
        ).execute()

    if (not args.freqtemp) and (options["solvent_model"] is not None or args.compound):
        # set or confirm solvent
        options["solvent"] = inquirer.fuzzy(
            message=f"Current solvent is {options['solvent']}. Press enter to confirm or type another solvent:",
            choices=epsilon_dict.keys(),
            default=options["solvent"],
            # validate=lambda result: result in epsilon_dict,
            # invalid_message="Solvent not recognized."
        ).execute()

def get_comp_script_inp(rootname):
    return f'''* xyzfile {options["charge"]} {options["mult"]} {rootname}.xyz
    
%pal
 nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

%geom
  MaxStep {options["maxstep"]}
end

%compound "{options["compound_job_scriptname"]}"
  with
    solvent = "{options["solvent"]}";
{options["comp_job_extra_variables"]}
end

'''

def get_freqtemp_inp(rootname):
    return f'''! PrintThermoChem

%geom
  inhessname "{rootname}.hess"
end

%freq Temp {options["temp"]} end

* xyzfile 0 1 {rootname}.xyz

'''

def get_solvent_block(solvent, solvent_model):

    if solvent_model == 'ALPB':
        return ""

    s = '%cpcm\n'
    if solvent_model == 'CPCM':
        s += f"  epsilon {epsilon_dict[solvent]}\n"
    elif solvent_model == 'SMD':
        s += f"  smd True\n  SMDsolvent \"{solvent}\"\n"
    else:
        raise Exception(f'Solvent model \"{solvent_model}\" not recognized (use CPCM, SMD or ALPB)')

    s += "end"

    return s

def get_inp(rootname, opt=True):
    maxstep_line = f'MaxStep {options["maxstep"]}'
    hess_line = "Calc_Hess true" if args.ts else ""

    if options["hh"]:
        hess_line += f'\n  Hybrid_Hess {"{"+hh_ids+"}"} end'

    geom_block = f'%geom\n\n  {maxstep_line}\n\n  {hess_line}\n\nend' if opt else ""

    solvent_keyword = {
        'ALPB' : f"ALPB({options['solvent']})",
        'CPCM' : 'CPCM',
        'SMD' : "",
    }[options['solvent_model']]

    inp = f'''! {options["level"]} {options["basis_set"]} {solvent_keyword} {options["opt"]}
! {options["additional_kw"]}

%pal
  nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

{geom_block}

{get_solvent_block(options["solvent"], options["solvent_model"])}

{f"%freq temp {options['temp']} end" if options["freq"] else ""}
{options["constr_block"]}
{options["extra_block"].replace('$ROOTNAME', rootname)}

* xyzfile {options["charge"]} {options["mult"]} {rootname}.xyz

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

    if args.goat:
        return 48

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

def inquire_constraints():

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

        indices = [int(i) for i in c_string.split()[:-1]]
        if c_type == "B" and c_string.split()[-1] == "ts":
            target = str(get_ts_d_estimate(args.inputfiles[0], indices))
            c_string = ' '.join(c_string.split()[:-1] + [target])
        else:
            target = float(c_string.split()[-1])

        constraint_strings.append((c_type, c_string)) 
        global_constraints.append(Choice(value=Constraint(indices, target), name=f"{c_type} {c_string}"))

    options["constr_block"] += "%geom Constraints\n"
    for c_type, c_string in constraint_strings:
        options["constr_block"] += "  {{ {0} {1} C }}\n".format(c_type, c_string)
        
    options["constr_block"] += "  end\nend"

def get_prev_popt_constr_ids(rootname):
    '''
    Looks in the parent folder for constraints in an ORCA .inp file.,
    and returns them as a list of InquirerPy Choices.

    '''
    choices = []

    popt_inp_name = rootname + '.inp'
    if popt_inp_name in os.listdir(os.path.dirname(os.getcwd())):
        popt_ids_lines = getoutput(f'grep {{*}} ../{popt_inp_name}')

        for line in popt_ids_lines.split("\n"):
            popt_ids = ' '.join(re.findall(r'(?<!\.)\b\d+\b(?!\.)', line))
            popt_ids = ' '.join({i for i in popt_ids.split()})

            line += ' - '
            for index in popt_ids.split():
                elem = getoutput(f'sed "{3+int(index)}q;d" {rootname}.xyz').split()[0]
                line += f'{elem}({index})-'
            line = line[:-1]
            choices.append(Choice(value=popt_ids, name=line))
        
        return choices

def inquire_ts_mode_following(default=""):

    letter_dict = {
        2 : "B",
        3 : "A",
        4 : "D",
    }

    choices = get_prev_popt_constr_ids(args.inputfiles[0][:-4])

    if choices:

        if inquirer.confirm(message=("Would you like the saddle opt. to follow the constrained "
                            f"coordinates from the previous popt?"), default=True).execute():

            popt_ids = inquirer.select(
                message="Select the constraint to follow:",
                choices=choices).execute()
            
            letter = letter_dict[len(popt_ids.split())]
            options["extra_block"] += "%geom\n  TS_Mode {{{0} {1}}}\n  end\nend".format(letter, popt_ids)
            return
    
    ask_manual = False
    if args.composite and global_constraints:
        ids_from_global = inquirer.select(
            message="Would you like the saddle opt. to follow one of these internal coordinates?",
                choices=global_constraints+[Choice(value="manual", name='Specify a different one'), Choice(value=False, name='None')],
                default=default,
            ).execute()

        if ids_from_global == 'manual':
            ask_manual = True

        elif not ids_from_global:
            return

        else:
            popt_ids = ' '.join([str(n) for n in ids_from_global.indices])
            options["extra_block"] += "%geom\n  TS_Mode {{{0} {1}}}\n  end\nend".format(ids_from_global.type, popt_ids)

    if ask_manual:
        manual_ids = inquirer.text(
            message=("Would you like the saddle opt. to follow a specific internal coordinate? "
                    f"If so, type the indices of the bond/angle/dihedral to follow."),
            default=default,
            ).execute()
    
    if manual_ids:
        popt_ids = manual_ids.rstrip().strip()
        letter = letter_dict[len(popt_ids.split())]
        options["extra_block"] += "%geom\n  TS_Mode {{{0} {1}}}\n  end\nend".format(letter, popt_ids)

def get_goat_block():

    s = "%goat\n"

    if options["level"] == "GFN2-XTB":
        s += "    GFNUPHILL GFNFF\n"

    s += "    MAXEN 6.0\n"
    s += "end\n\n"

    return s

def inquire_level(default='wB97M-V'):

    options["level"] = inquirer.select(
            message="Which level of theory would you like to use?:",
            choices=(
                Choice(value='B3LYP D3BJ', name='B3LYP-D3(BJ) - GGA           ★★★☆☆'),
                Choice(value='PBE0',       name='PBE0         - GGA           ★★★☆☆'),
                Choice(value='R2SCAN-3c',  name='R²SCAN-3c    - (composite)   ★★★☆☆'),
                Choice(value='wB97X-V',    name='ωB97X-V      - GGA           ★★★★☆'),
                Choice(value='wB97X-D4',   name='ωB97X-D4     - GGA           ★★★★☆'),
                Choice(value='wB97M-V',    name='ωB97M-V      - Meta-GGA      ★★★★☆'),
                Choice(value='M062X',      name='M06-2X       - Meta-GGA      ★★★★☆'),
                Choice(value='wB97M(2)',   name='ωB97M(2)     - Double Hybrid ★★★★★ - (will copy .gbw from parent folder if not in current)'),
            ),
            default=default,
        ).execute()

def set_auto_solvent():
    '''
    Modifies the options["solvent"] value to the most popular
    if input files are found the current folder.
    
    '''

    solvent_list = []

    # look for CPCM solvents via epsilon
    try:
 
        epsilons = [float(line.split()[1]) for line in getoutput('grep epsilon *.inp -h').split('\n')]
        for e in epsilons:
            solvent_list.append(
                next(key for key, value in epsilon_dict.items() if value == e)
            )
    except (IndexError, ValueError, StopIteration):
        pass

    # look for SMD solvents
    try:
        for line in getoutput('grep SMDsolvent *.inp -h').split('\n'):
            solvent = line.split()[1].lstrip('\"').rstrip('\"')
            if solvent in epsilon_dict.keys():
                solvent_list.append(solvent)
        
    except (IndexError, ValueError, StopIteration):
        pass

    # look for solvents in compound jobs
    try:
        for line in getoutput('grep solvent *.inp -h').split('\n'):
            solvent = line.split()[2].rstrip(";").lstrip('\"').rstrip('\"')
            if solvent in epsilon_dict.keys():
                solvent_list.append(solvent)
        
    except (IndexError, ValueError, StopIteration):
        pass

    if not solvent_list:
        return

    options["solvent"] = sorted([(n, solvent_list.count(n)) for n in set(solvent_list)], key=lambda x: x[1], reverse=True)[0][0]

def copy_from_parent_if_not_here(filename):
    '''
    Copy a file from the parent folder to the current,
    unless a file with the same name is already present
    in the current folder.
    
    '''

    if not filename in os.listdir():
        parent = os.path.dirname(os.getcwd())
        if filename in os.listdir(parent):
            source = os.path.join(parent, filename)
            target = os.path.join(os.getcwd(), filename)
            shutil.copyfile(source, target)
            print(f'--> Copied {filename} from parent folder to current folder.')
        else:
            raise Exception(f'Could not find {filename} file in current nor parent folder.')
        
def multiplicity_check(rootname, charge, multiplicity=1) -> bool:
    '''
    Returns True if the multiplicity and the nuber of
    electrons are one odd and one even, and vice versa.

    '''

    electrons = 0
    for line in getoutput(f'cat {rootname}.xyz').splitlines():
        parts = line.split()
        if len(parts) == 4:
            try:
                element = parts[0]
                electrons += getattr(pt, element).number
            except AttributeError:
                pass

    electrons -= charge
    
    return (multiplicity % 2) != (electrons % 2)
        
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
parser.add_argument("-priority", help="Run jobs with priority (and estimate cost).", action="store_true", required=False)
parser.add_argument("-goat", help="Conformational search via GOAT.", action="store_true", required=False)
parser.add_argument("-freqtemp", help="Recalculate vibrational corrections at a new temperature.", action="store_true", required=False)
args = parser.parse_args()

if not args.inputfiles:
    raise Exception('No input structures selected. Please specify at least one.')

set_auto_solvent()

if not any((args.sp,
            args.fastsp,
            args.ts,
            args.optf,
            args.popt,
            args.scan,
            args.nmr,
            args.compound,
            args.irc,
            args.goat)):
    inquirer_set_options(args)

if args.solvent:
    print(f"--> Setting solvent to {options['solvent']}")

if args.sp:

    inquire_level()

    options["opt"] = ""
    options["freq"] = False
    options["basis_set"] = 'def2-TZVPP'
    # options["solvent_model"] = 'SMD'
    options["additional_kw"] += " Defgrid3"

    if options["level"] == 'wB97M(2)':
        options["additional_kw"] += " SCNL NoFrozenCore def2-TZVPP/C CALCGUESSENERGY NOITER"
        options["extra_block"] += "%scf\n  Guess MORead\n  MOInp \"$ROOTNAME.gbw\"\nend"

if args.fastsp:
    options["opt"] = ""
    options["freq"] = False
    options["level"] = 'R2SCAN-3c'
    options["basis_set"] = ''
    options["additional_kw"] += " Defgrid3"

if args.ts:
    options["opt"] = "OptTS"
    options["mem"] = 8
    options["freq"] = True
    options["additional_kw"] += " TightOpt LARGEPRINT Defgrid3"

    inquire_ts_mode_following()

if args.hh:
   
    hh_ids = input("Hybrid Hessian: provide indices involved in TS: ")

    if hh_ids == "":

        popt_ids = get_prev_popt_constr_ids(args.inputfiles[0][:-4])
        
        if popt_ids:
            if inquirer.confirm(f"Would you like to use indices from the previous popt? [{popt_ids}]", default=True):
                hh_ids = popt_ids

            else:
                print('Please specify hybrid hessian-excluded indices.')
                sys.exit()
        else:
            print('Please specify hybrid hessian-excluded indices.')
            sys.exit()

if args.optf:
    options["opt"] = "Opt"

    options["freq"] = inquirer.confirm(message='Run frequency calculation?', default=False).execute()
    if inquirer.confirm(message='LARGEPRINT? (writes orbital information to .out)', default=False).execute():
        options["additional_kw"] += " LARGEPRINT"
    options["additional_kw"] += " TightOpt Defgrid3"

if args.popt:
    options["freq"] = False
    options["opt"] = "Opt"
    options["popt"] = True
    options["additional_kw"] += " Defgrid3"

    inquire_constraints()

if args.scan:
    options["opt"] = "Opt"

    # make scan %geom block and add it to the extra block
    options["constr_block"] += get_scan_block(args.inputfiles[0])

if args.nmr:
    options["freq"] = False
    options["opt"] = ""
    options["level"] = 'PBE0'
    options["basis_set"] = '6-311+G(2d,p)'
    options["additional_kw"] += " NMR Defgrid3"

if args.compound:

    options["opt"] = True
    options["compound_job_scriptname"] = inquirer.select(
        message='What compound script would you like to run?',
        choices=(
            Choice(value="optf+sp.cmp", name='optf+sp.cmp - Unconstrained opt., freq. calc., single point energy calc.'),
            Choice(value="popt+saddle+sp.cmp", name='popt+saddle+sp.cmp  - Constrained opt., saddle point opt., freq. calc., single point energy calc.'),
        ),
        default="optf+sp.cmp",
    ).execute()

    try:
        script_folder = os.path.dirname(os.path.realpath(__file__))
        source = os.path.join(script_folder, options["compound_job_scriptname"])
        shutil.copyfile(source, options["compound_job_scriptname"])
        print(f'--> Copied {options["compound_job_scriptname"]} to input file folder.')
    except Exception:
        raise Exception(f'Something went wrong when copying {options["compound_job_scriptname"]} to the destination folder.')
    
    if options["compound_job_scriptname"] == "popt+saddle+sp.cmp":
        inquire_constraints()
        inquire_ts_mode_following()
        options["comp_job_extra_variables"] = f'    extra_block_1 = \"{options["constr_block"]}\";\n    extra_block_2 = \"{options["extra_block"]}\";\n'
        options["freq"] = True

if args.irc:
    options["freq"] = True
    options["additional_kw"] += " IRC Defgrid3"
    options["opt"] = ""

if args.goat:
    options["level"] = inquirer.select(
            message="What level of theory would you like to use?",
            choices=(
                Choice(value='GFN-FF', name='GFN-FF'),
                Choice(value='GFN2-XTB', name='GFN2-XTB (will still request GFN-FF for the uphill steps)'),
            ),
            default='GFN-FF',
        ).execute()

    options["opt"] = ""
    options["additional_kw"] += " GOAT"
    options['solvent_model'] = 'ALPB'
    options["extra_block"] += get_goat_block()

    if inquirer.confirm(
                message="Do you wish to specify any constraint?",
                default=False,
            ).execute():
        inquire_constraints()

if args.freqtemp:
    ans = inquirer.text(message="Specify the new temperature for the vibrational correction, in degrees Celsius:").execute()
    options["temp"] = round(float(ans) + 273.15, 2)

# Make sure to add Freq for appropriate calculations if we have not already
if options["freq"] and "Freq" not in options["additional_kw"]:
    options["additional_kw"] += " Freq"

# Inquire about automatic charges based on file names
allchars = ''.join(args.inputfiles)
auto_charges = False
if '+' in allchars or '-' in allchars:
    auto_charges = inquirer.confirm(
        message="Found charge signs in filenames (+/-). Auto assign charges in input files?",
        default=False,
    ).execute()

# print options
print()
for option, value in options.items():
   if value != "":
       print(f"--> {option} = {value}")
print()

# calculate processors to use for each job and ask if user wants to override
if args.freqtemp:
    procs_list = [1 for filename in args.inputfiles]

else:
    procs_list = [get_procs(filename.split('.')[0]) for filename in args.inputfiles]
    if newprocs := inquirer.text(
            message=f'Automatic number of cores is {procs_list[0]}. Type a different number to override or enter to confirm.'
                ).execute():
        procs_list = [int(newprocs) for _ in procs_list]

cum_cost, cum_cpu, cum_mem = 0, 0, 0

#################################################################################### COPY STUFF FROM PARENT FOLDER

# make sure we have all the required files for the job in the submission folder

# wB97M(2) needs the gbw file (copy from parent folder)
if options["level"] == 'wB97M(2)':
    for filename in args.inputfiles:
        wfn = filename[:-4]+".gbw"
        copy_from_parent_if_not_here(wfn)
        
# freqtemp needs the hessian (copy from parent folder)
if args.freqtemp:
    for filename in args.inputfiles:
        hessname = filename[:-4]+".hess"
        copy_from_parent_if_not_here(hessname)

#################################################################################### WRITE INPs, PRINT JOB TABLE

# start printing job table
print(f'#    Filename                       Cores  Mem(GB)  Max cost (24h)')
print( '------------------------------------------------------------------')

for f, (filename, procs) in enumerate(zip(args.inputfiles, procs_list)):
    rootname = filename.split('.')[0]
    options["procs"] = procs
    cost = get_daily_cost(procs, mem_per_core=options["mem"])
    cum_cost += cost
    cum_cpu += procs
    cum_mem += options["mem"] * procs

    print(f'{str(f+1):>3}  {filename:30s} {procs:2}     {options["mem"]*procs:d}       {cost:.2f} $')

    if auto_charges:
        if '+' in rootname:
            options['charge'] = 1
            print(f'--> {filename} : assigned charge +1 based on input name')
        elif '-' in rootname:
            options['charge'] = -1
            print(f'--> {filename} : assigned charge -1 based on input name')
        else:
            options['charge'] = 0
            print(f'--> {filename} : assigned charge 0 based on input name')

    if multiplicity_check(rootname, options["charge"]):
        options['mult'] = 1
    else:
        options['mult'] = inquirer.text(
            message=f'It appears {rootname} is not a singlet. Please specify multiplicity:',
            validate=lambda inp: inp.isdigit() and int(inp) > 1,
        ).execute()

    if args.compound:
        s = get_comp_script_inp(rootname)

    elif args.freqtemp:
        s = get_freqtemp_inp(rootname)

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

print( '------------------------------------------------------------------')
print(f'Maximum estimated cost (24 h runtime, {int(cum_cpu)} CPUs, {int(cum_mem)} GB MEM): {cum_cost:.2f} $\n')

run_jobs = inquirer.confirm(
        message="Run jobs for the newly generated input files?",
        default=True,
    ).execute()

if run_jobs:
    inpnames = [name.split(".")[0]+".inp" for name in args.inputfiles]

    from orcasub_batch import main

    # adding dummy first element to simulate sys.argv
    main([None]+inpnames)