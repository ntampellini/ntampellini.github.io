---
layout: post
title: Scripts
---
<style> .myButton {
  background-color:#ffd9ce;
	border-radius:17px;
	display:inline-block;
	cursor: pointer;
	color:#000000;
	/* font-family: Arial; */
	/* font-size: 28px; */
	padding: 7px 12px;
	text-decoration: none;
}
.myButton:hover {
    background-color:#e86f87;
}
.myButton:active {
    position: relative;
	top:1px;
} </style>

Here are some Python and Bash scripts that I routinely use and share, divided by use case. In most cases, running these "dry" (with no arguments) prints the usage text you can see here. Click on the button to download a script.

# ORCA/XTB

<a href="/assets/scripts/utils.py" class="myButton"
   title="Download script file">utils.py</a>

Small library of Python functions required by some of these scripts.

<br>

<a href="/assets/scripts/orcasub.sh" class="myButton"
   title="Download script file">orcasub.sh</a>

    Slurm batch script to run ORCA calculations with the latest setup for
    Grace, Yale's HPC cluster. Usage: 

      orcasub.sh input[.inp/.xyz] [n_cores]

    If ncores is not specified, it is tentatively extracted from input.inp

<br>

<a href="/assets/scripts/mkorca.py" class="myButton"
   title="Download script file">mkorca.py</a>

    Makes one or more ORCA inputs following the desired specifications. Syntax:

      python mkorca.py conf*.xyz [option]

    conf*.xyz: base name of input geometry file(s)
    option:
      sp:     single-point energy calculation
      fastsp: fast single-point energy calculation
      optf:   optimization + frequency calculation
      popt:   partial optimization (specify distance constraint)
      ts:     eigenvector-following saddle point optimization
      nmr:    single-point energy calculation with NMR shieldings.

    Each of these options might have different default levels of theory.
    Manually check/modify this script at your convenience.

<br>

<a href="/assets/scripts/mkorcascan.py" class="myButton"
   title="Download script file">mkorcascan.py</a>

    Makes an ORCA input file for a distance scan. Syntax:

      python mkorcascan.py filename.xyz i1 i2 target

    filename.xyz: base name of input geometry file
    i1/i2: indices of atoms to scan
    target: desired distance between indices at end of scan

<br>

<a href="/assets/scripts/orcasub_batch.py" class="myButton"
   title="Download script file">orcasub_batch.py</a>

    Launch batches of ORCA jobs through orcasub.sh. Syntax:

      python orcasub_batch.py conf*.xyz [n]

    n: optional, overrides the default value of 16 for cores

<br>

<a href="/assets/scripts/check.py" class="myButton"
   title="Download script file">check.py</a>

    Check running and completed ORCA jobs and prints/plots relevant information.
    For scans, also extracts local maxima structures. Syntax:

      python check.py filename[.out] [i1] [i2] [zoom]

    filename: base name of input/output file (with or without extension)
    i1/i2: optional, plot i1-i2 distance during optimization steps
    zoom: optional, only shows plot of last 20 iterations

<br>

<a href="/assets/scripts/compare.py" class="myButton"
   title="Download script file">compare.py</a>

    Compare completed ORCA jobs electronic/free energy of groups of jobs and
    prints collected energy values. Syntax:

      python compare.py conf*.out [g] 

    conf*.out: base name of output file(s)
    g: optional, uses free energy for the comparison (needs freq calculations)

<br>

<a href="/assets/scripts/update.py" class="myButton"
   title="Download script file">update.py</a>

    Updates ORCA input files by replacing them with the last step of
    the optimization trajectory. Syntax:

      python update.py conf*.xyz

    conf*.xyz: base name of input geometry file(s)
       
<br>

<a href="/assets/scripts/xtbopt.py" class="myButton"
   title="Download script file">xtbopt.py</a>

    Optimizes the specified geometry/ies, compares results and replaces
    the input file. Syntax:

      python xtbopt.py filename*.xyz [newfile] [ff] [sp] [c] [charge=n]
    
    filename*.xyz: Base name of input geometry file(s)
    newfile: Optional, creates a new file, preserving the original
    ff: Optional, use GFN-FF instead of GFN2-XTB
    sp: Optional, run a single point energy calculation
    c: Optional, specify one or more distance/dihedral constraints
    charge=n: Optional, where "n" is an integer. Specify the total charge

<br>