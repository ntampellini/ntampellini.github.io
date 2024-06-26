####################################################
options = {
    "solvent" : "toluene",
    "solvent_model" : 'CPCM',
    "level" : "R2SCAN-3c",
    "basis_set" : "",
    "opt": "",
    "ts" : False, # Set true for saddle optimization
    "freq" : False,
    "temp" : 273.15+25,
    "procs" : 16,
    "mem" : 6, # Memory per core, in GB
    "charge" : 0,
    "maxstep" : 0.05, # in Bohr atomic units (1au = 0.529177 A)
    "popt" : False,
    "compound" : False,

    "additional_kw" : "Defgrid3",
    "extra_block" : "",
    }

epsilon = {
     "dichloromethane" : 9.04,
     "chloroform" : 4.81,
     "phcf3" : 9.18,
     "dmso" : 47.2,
     "toluene" : 2.4,
}

# round temperature so it looks prettier
options["temp"] = round(options["temp"], 2)

####################################################

import os 
import sys

if len(sys.argv) == 1:
    print("\n  Makes one or more ORCA inputs following the desired specifications. Syntax:\n\n" +
          "  python mkorca.py conf*.xyz [option]\n\n" + 
          "  conf*.xyz: base name of input geometry file(s)\n" +
          "  option:\n" +
          "    sp:     single-point energy calculation.\n" +
          "    fastsp: fast single-point energy calculation.\n" +
          "    optf:   optimization + frequency calculation.\n" +
          "    popt:   partial optimization (specify distance constraint).\n" +
          "    ts:     eigenvector-following saddle point optimization.\n" +
          "    nmr:    single-point energy calculation with NMR shieldings.\n" + 
          "  Each of these options might have different default levels of theory. Manually check/modify this script at your convenience.\n"
        )
    quit()

def get_comp_script_inp(rootname):
    return f'''* xyzfile 0 1 {rootname}.xyz
    
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
        s += f"  epsilon {epsilon[options['solvent']]}\n"
    elif solvent_model == 'SMD':
        s += f"  smd True\n  SMDsolvent \"{options['solvent']}\"\n"
    else:
        raise Exception(f'Solvent model \"{solvent_model}\" not recognized (use CPCM or SMD)')

    s += "end"

    return s

def get_inp(rootname):
    return f'''! {options["level"]} {options["basis_set"]} {"CPCM" if options["solvent_model"] == "CPCM" else ""} {options["opt"]}
! {options["additional_kw"]}

%pal
  nprocs {options["procs"]}
end

%maxcore {options["mem"]*1024}

%geom
  MaxStep {options["maxstep"]}
  {"Calc_Hess true" if options["ts"] else ""}
end

{get_cpcm_block(options["solvent"], options["solvent_model"])}

{f"%freq temp {options['temp']} end" if options["freq"] else ""}

{options["extra_block"]}

* xyzfile {options["charge"]} 1 {rootname}.xyz

'''

os.chdir(os.getcwd())
names = []

if "sp" in sys.argv:
    sys.argv.remove("sp")
    options["ts"] = False
    options["opt"] = ""
    options["procs"] = 16
    options["freq"] = False
    options["level"] = 'wB97M-V'
    options["basis_set"] = 'def2-TZVPP'
    options["solvent_model"] = 'SMD'

if "fastsp" in sys.argv:
    sys.argv.remove("fastsp")
    options["ts"] = False
    options["opt"] = ""
    options["procs"] = 16
    options["freq"] = False
    options["level"] = 'B97-3c'
    options["basis_set"] = ''

if "ts" in sys.argv:
    sys.argv.remove("ts")
    options["ts"] = True
    options["opt"] = "OptTS"
    options["procs"] = 16
    options["freq"] = True
    options["additional_kw"] += " TightOpt LARGEPRINT"

if "optf" in sys.argv:
    sys.argv.remove("optf")
    options["ts"] = False
    options["opt"] = "Opt"
    options["procs"] = 16
    options["freq"] = True
    options["additional_kw"] += " TightOpt LARGEPRINT"

if "popt" in sys.argv:
    sys.argv.remove("popt")
    options["ts"] = False
    options["freq"] = False
    options["opt"] = "Opt"
    options["popt"] = True
    options["procs"] = 16
    options["additional_kw"] += " TightOpt LARGEPRINT"
    c_string = input("Provide indices to constrain and distance: ")
    assert len(c_string.split()) == 3

    options["extra_block"] = "%geom Constraints\n  {{ B {0} C }}\n  end\nend".format(c_string)

if "nmr" in sys.argv:
    sys.argv.remove("nmr")
    options["ts"] = False
    options["freq"] = False
    options["opt"] = ""
    options["procs"] = 16
    options["level"] = 'PBE0'
    options["basis_set"] = '6-311+G(2d,p)'
    options["additional_kw"] += " NMR"

if "comp" in sys.argv:
    options["compound"] = True
    options["procs"] = 8

if "irc" in sys.argv:
    sys.argv.remove("irc")
    options["irc"] = True
    options["ts"] = False
    options["freq"] = True
    options["procs"] = 32
    options["additional_kw"] += " IRC"
    options["opt"] = ""

if options["freq"] and"Freq" not in options["additional_kw"]:
    options["additional_kw"] += " Freq"


for option, value in options.items():
   if value != "":
       print(f"--> {option} = {value}")


for filename in sys.argv[1:]:
    rootname = filename.split('.')[0]

    if options["compound"]:
        s = get_comp_script_inp(rootname)

    else:
        s = get_inp(rootname)

        # if options["irc"]:
        #     options["extra_block"] = "%irc\n  InitHess read\n  Hess_Filename \"{0}.hess\"\nend\n".format(rootname)

    with open(f'{rootname}.inp', 'w') as f:
        f.write(s)

    print(f'Written orca input file {rootname}.inp')

    # Convert all text files to Linux format
    os.system(f'dos2unix {rootname}.inp')
    os.system(f'dos2unix {rootname}.xyz')

print('Converted text files to Linux format.')