---
layout: post
title: Energy Converter
---

<HTML>
<HEAD>
<TITLE>Energy Converter</TITLE>

<p>Energy converter tool adapted from <a href="https://www.colby.edu/chemistry/PChem/Hartree.html">Colby College</a>.</p>

<SCRIPT LANGUAGE="JavaScript">

//DEFINE METHODS
function constants(conv) {
    var numE = 7;
          conv[0] = 1.000 ;
// eV
          conv[1] = 2.7211399E+01 ;
// kJ/mol
          conv[2] = 2.6255002E+03 ;
// kcal/mol
          conv[3] = conv[2]/4.184 ;
// cm-1
          conv[4] = 2.1947463E+05 ;
// V
          conv[5] = 2.6255002E+06/96484.6 ;
// K
          conv[6] =3.1577709E+05 ;
// Boltzmann
          conv[7] = -conv[6] ;
            return numE;
}

function springConstants(conv) {
   var numE = 4;
// Eh/bohr**2
          conv[0] = 1.000 ;
// Eh/A**2
          conv[1] = 1/((0.529177249)**2) ;
// eV/A**2
          conv[2] = 2.7211399E+01 * conv[1] ;
// kcal/A**2
          conv[3] = 627.5096080305927 * conv[1] ;
   return numE;
}

function displayInfo(form,field) {
// find field index

   for (var i=0; i<=(spring_nfields+nfields+1); i++) {
      if ( form.elements[i].name == field ) {
         idx = i ;
         break;
      }
   }
   // console.log("idx is", idx)
   if (idx >= 9) { displaySpringInfo(form,field) } else {

      // find number of characters in input string for significant figure functions
      nchars = form.elements[idx].value.length +1 ;

      // calculate the base energy in Hartrees
      if ( idx != 7) {
         energy = form.elements[idx].value/conv[idx];
      } else {
         energy = Math.log(form.elements[idx].value)*298.15/conv[idx];
      }

      // convert to other units
      for (var i=0; i<=nfields; i++) {
         if ( i != idx ) {
            if ( i != 7) {
               form.elements[i].value = trunc(energy*conv[i],nchars) ;
            } else {
               form.elements[i].value = trunc(Math.exp(energy*conv[i]/298.15),4) ;
            }
         }
      }

      boltzmann()
   }
}

function setT() {
    var num = Number(document.Boltzmann.C.value) + 273.15;
    nchars = Math.floor(num).toString().length + 3; 
    document.Boltzmann.T.value = trunc(num, nchars); 
}

function setC() {
    var num = Number(document.Boltzmann.T.value) - 273.15;
    nchars = Math.floor(num).toString().length + 3;
    document.Boltzmann.C.value = trunc(num, nchars); 
}

function boltzmann() {
// calculate boltzmann fractions and voltage for general conditions
      var T = document.Boltzmann.T.value ;
      // var C = document.Boltzmann.C.value ;
      var gj = document.Boltzmann.gj.value ;
      var gi = document.Boltzmann.gi.value ;
      var z= document.Boltzmann.z.value ; 
      var r = Math.exp(energy*conv[7]/T)*gj/gi ;
      document.Boltzmann.flow.value = trunc(1/(r+1)*100.0,3) ;
         var fup = trunc(r/(r+1)*100.0,3) ;
    if ( fup > 1e-20 ) {
             document.Boltzmann.fup.value = fup 
                              } else {
             document.Boltzmann.fup.value = 0 }
// Put in diagram poulations
         //  var all = "-oooooooooo" ;
         //  var molecules = Math.floor(r/(r+1)*10.0+0.5) ;
         //    document.Boltzmann.up.value = all.substring(0,molecules+1) + "-" ;
         //    document.Boltzmann.low.value = all.substring(0,11-molecules) + "-" ;
// Voltage for z != 1
      voltage = energy*conv[5]/z
      nchars = Math.floor(voltage).toString().length + 3;
      document.Boltzmann.V.value = trunc(voltage, nchars);

      halflife();
}

function halflife() {
   var T = document.Boltzmann.T.value;
   var gj = document.Boltzmann.gj.value;
   var gi = document.Boltzmann.gi.value;
   var kb = 1.380649E-23;
   var h = 6.62607015E-34;
   var R = 0.001985877534;
   k_rate = kb * T / h * Math.exp(-energy*conv[3]/(R*T));
   time = Math.log(2) / k_rate;
   
   if (time > 3600*24*365) {
      uom = "years";
      factor = 3600*24*365;
   } else if (time > 3600*24) {
      uom = "days";
      factor = 3600*24;
   } else if (time > 3600) {
      uom = "hours";
      factor = 3600;
   } else if (time > 60) {
      uom = "minutes";
      factor = 60;
   } else if (time > 1) {
      uom = "seconds";
      factor = 1;
   } else if (time > 1E-3) {
      uom = "ms";
      factor = 1E-3;
   } else if (time > 1E-6) {
      uom = "μs";
      factor = 1E-6;
   } else if (time > 1E-9) {
      uom = "ns";
      factor = 1E-9;
   } else {
      uom = "ps";
      factor = 1E-12;
   }

   console.log(time/factor, uom)
   intlength = Math.floor(time/factor).toString().length
   if (intlength > 8) {
      nchars = 4;
   } else {
      nchars = intlength + 3;
   }
   document.Boltzmann.halflife.value = trunc(time/factor, nchars).toString() + " " + uom;
}

// Significant figure functions
function ord(x) {
   return Math.floor(Math.log(Math.abs(x+1e-35))/2.303)
}
// Truncate to n sign. figures
function trunc(x,n) {
   return Math.floor(x*Math.pow(10,-ord(x)+n-1)+.5)/Math.pow(10,-ord(x)+n-1)
}

function displaySpringInfo(form,field) {

   // find field index
   for (var i=0; i<=(spring_nfields+nfields+1); i++) {
      if ( form.elements[i].name == field ) {
         idx = i ;
         break;
     }
   }

   // find number of characters in input string for significant figure functions
   nchars = form.elements[idx].value.length +3 ;

   // calculate the base spring force constant in Hartrees over A squared
   k = form.elements[idx].value/springConv[idx-9];
   console.log("idx is", idx)
   console.log("k is", k)

   // convert to other units
   for (var i=0; i<=(spring_nfields+nfields+1); i++) {
      if ( i != idx ) {
         form.elements[i].value = trunc(k*springConv[i-9],nchars) ;
      }
      }
   }

// MAIN variable declarations
 var energy = 0.000;
 var nchars = 0;
 var conv = new Array();
 var nfields = constants(conv);

 var k = 0.000;
 var springConv = new Array();
 var spring_nfields = springConstants(springConv)

</SCRIPT>
</HEAD>
<BODY>
<!-- <H3>Energy Units Converter</H3> -->
Enter your energy value in the box with the appropriate units, then press "tab"
or click outside of the input box.
<P>
<FORM NAME="Hartree" METHOD="POST">
<INPUT TYPE="text" NAME="H" VALUE="0" onChange="displayInfo(this.form,this.name);"> Hartrees (Eh)<BR>
<INPUT TYPE="text" NAME="eV" VALUE="0" onChange="displayInfo(this.form,this.name);"> eV<BR>
<INPUT TYPE="text" NAME="kJ/mol" VALUE="0" onChange="displayInfo(this.form,this.name);"> kJ/mol<BR>
<INPUT TYPE="text" NAME="kcal/mol" VALUE="0" onChange="displayInfo(this.form,this.name);"> kcal/mol<BR>
<INPUT TYPE="text" NAME="cm-1" VALUE="0" onChange="displayInfo(this.form,this.name);"> cm<sup>-1</sup><BR>
<INPUT TYPE="text" NAME="V" VALUE="0" onChange="displayInfo(this.form,this.name);"> V for 1e<sup>-</sup> transfer<BR>
<INPUT TYPE="text" NAME="K" VALUE="0" onChange="displayInfo(this.form,this.name);"> K (equivalent temperature)<BR>
<INPUT TYPE="text" NAME="B" VALUE="1" onChange="displayInfo(this.form,this.name);"> Boltzmann population ratio at 25 °C
</FORM>
<HR>
<H3>Boltzmann Populations and kinetics</H3>
Change the temperature and degeneracies:
<FORM NAME="Boltzmann" METHOD="POST">
<TABLE>
<TR><TD><INPUT TYPE="text" NAME="T" SIZE="5" VALUE="298.15" onChange="setC();boltzmann();"> K / <INPUT TYPE="text" NAME="C" SIZE="5" VALUE="25" onChange="setT();boltzmann();"> °C</TD>
<TD></TD>
<!-- <TD></TD> -->
</TR>
<TR><TD>
upper degeneracy: <INPUT TYPE="text" NAME="gj" SIZE="2" VALUE="1" onChange="boltzmann();"></TD>
<TD><INPUT TYPE="text" NAME="fup" SIZE="7" VALUE="50.0" onFocus="blur();"> %</TD>
<!-- <TD><INPUT TYPE="text" NAME="up" SIZE="12" VALUE="-ooooo-" onFocus="blur();"></TD> -->
</TR>
<TR><TD>
lower degeneracy: <INPUT TYPE="text" NAME="gi" SIZE="2" VALUE="1" onChange="boltzmann();"></TD>
<TD><INPUT TYPE="text" NAME="flow" SIZE="7" VALUE="50.0" onFocus="blur();"> %</TD>
<!-- <TD><INPUT TYPE="text" NAME="low" SIZE="12" VALUE="-ooooo-" onFocus="blur();"></TD> -->
</TR>
</TABLE>

Half-life time (first-order kinetic regime): <INPUT TYPE="text" NAME="halflife" VALUE="0" onFocus="blur();">

<H3>Electrochemical Voltage</H3>
Change the number of electrons transferred, z=
<INPUT TYPE="text" NAME="z" SIZE="2" VALUE="1" onChange="boltzmann();"><BR>
Then, V=
<INPUT TYPE="text" NAME="V" SIZE="7" VALUE="0" onFocus="blur();"> V<BR>

<H3>Spring constants</H3>
<INPUT TYPE="text" NAME="Eh_bohr" VALUE="0" onChange="displayInfo(this.form,this.name);"> Eh/bohr<sup>2</sup><BR>
<INPUT TYPE="text" NAME="Eh_A" VALUE="0" onChange="displayInfo(this.form,this.name);"> Eh/A<sup>2</sup><BR>
<INPUT TYPE="text" NAME="eV_A" VALUE="0" onChange="displayInfo(this.form,this.name);"> eV/A<sup>2</sup><BR>
<INPUT TYPE="text" NAME="kcal_A" VALUE="0" onChange="displayInfo(this.form,this.name);"> kcal/A<sup>2</sup><BR>


</FORM>