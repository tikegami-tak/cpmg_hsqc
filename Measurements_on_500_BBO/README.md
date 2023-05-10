Each folder saves a pulse programs (with the name beginning with ti) and parameter files (Topspin 3.6.2). The four folders correspond to the four types of CPMG methods. These are the data measured by Burker 500 MHz NMR (with a room-temperature BBO probe), but the parameters for a cryo-probe are almost the same except for 13C and 1H RF powers.

The description of the experiment can be found in pdata/1/title.

First, assemble the pulse program part of the CPMG period using the ti13Ccpmg_hsqc4.pl program. Then, paste the output into the specified location in the pulse program ti13Ccpmg_hsqc*. The Perl and pulse programs are saved in each folder.

After the measurement, separate the "ser" file with another Perl program "sepser_N.pl". In the program, specify the number of points in the direct measurement dimension (e.g., 2,048) and the number of files (e.g., 18, equal to NBL) to split into. The divided files have names like ser_1, ser_2, ..., ser_18.

Then go down to the folder "process" and run the Perl program "all_convft1.pl". It tells the "convft.com" program to read the separated ser_1, ser_2, etc. in order and to perform the Fourier transform.

To assign the peak(s), you can use any software, such as Sparky.

We used nmrDraw to automatically detect peaks from the reference spectrum only. The program "draw_peak3_noassign.pl" in the folder "pick_area" reads out the coordinates of, for example, 14 pixels centered on the peak top. The coordinates are saved in the file "extpoint.dat". The positions of the 14 pixels can be confirmed by displaying the file "temp.ft" through nmrDraw.

If everything looks good, run "intense_plot4.pl", which extracts the intensities at the same 14 pixel positions from the other spectra and saves them into a file. The intensity data are saved in the folder "intense_plot" and can be visualized with a command "gnuplot 1X.plot".

Similarly,  run "reff_plot6.pl", which calculates the effective R2 and saves them into a file. The effective R2 data are saved in the folder "reff_plot" and can be visualized with a command "gnuplot 1X.plot".
