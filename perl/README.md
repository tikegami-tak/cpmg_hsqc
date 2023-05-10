The Perl and Mathematica programs solved the Bloch-McConnell equation, whereas this Perl program solves the density matrix equation and includes no exchange phenomenon.

Furthermore, it currently implements only the conventional RC and AFTAC methods, which can be switched by setting the variable "sem" to 3 or 2, respectively.

The path to the library must be set correctly and written on the first line of the program (e.g., -I../TI_lib_perl). Then you can simply type "cpmg_hsqc_CH3_6.pl" in the shell to run it.

The graph can be seen by typing "gnuplot gnu_sem2.plot".
