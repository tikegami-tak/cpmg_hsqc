The content of this Mathematica program is basically the same as that of the Python program. However, the Mathematica program has no fitting function.

You can use the paid Mathematica or the free-of-charge Wolfram Engine running on Jupyter lab.

Change the variables written from lines 141 to 215 to appropriate values. The CPMG method can be changed with "hpi" on line 10. The meanings of the variables are the same as those in the Python program.

Since variables in Mathematica programs are global by default, you need to erase all the previous variables by typing "Remove["Global`@*"]" every time you run a program. After that, load the program with "<< cpmg_hsqc_exch_21.m".

You can draw a CPMG curve by running the "bloch[3]" function.
