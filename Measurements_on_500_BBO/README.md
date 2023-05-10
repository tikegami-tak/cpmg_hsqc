Each folder saves a pulse programs (beginning with ti) and parameter files (Topspin 3.6.2). The 4 folders correspond to the 4 types of CPMG methods, respectively. These are the data when measured by Burker 500 MHz NMR (BBO probe), but the parameters for a cryo-probe are almost the same except for 13C and 1H RF powers.

First, assemble the pulse program portion of the CPMG period using the ti13Ccpmg_hsqc4.pl program. Then, paste the output into the specified location in the pulse program ti13Ccpmg_hsqc*. The Perl and pulse programs are saved in each folder.

After the measurement, separate the "ser" file with another Perl program "sepser_N.pl". In the program specify the number of points in the direct measurement dimension (e.g., 2,048) and the number of files (e.g., 18) to split into. The divided files have names like ser_1, ser_2, ..., ser_18.

Then go to the folder "process" and run the Perl program "all_convft1.pl". It makes the "convft.com" program read the separated ser_1, ser_2, etc. in order and perform the Fourier transform.

To assign the peak(s), you can use any software, such as Sparky.

We used nmrDraw to automatically detect peaks from the reference spectrum only. The program "draw_peak3_noassign.pl" in the folder "pick_area" reads out the coordinates of, for example, 14 pixels centered on the peak top. The coordinates are saved in the file "extpoint.dat". The positions of the 14 pixels can be confirmed by displaying the file "temp.ft" through nmrDraw.

If everything looks good, run "intense_plot4.pl", which extracts the intensities at the same 14 pixel positions from the other spectra and saves them into a file. The file is saved in the folder "intense_plot" and can be visualized with "gnuplot 1X.plot".

I'll continue soon.
