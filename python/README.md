The Python program simulates the four types of relaxation dispersion curves, as mentioned in the upper folder.

The following assumes that you are using jupyter lab.

To execute the program, type “run cpmg_hsqc_exch_fit_34.py” and type one of the functions listed at the bottom of the program.
When you have modified the program, type “reset -f” before running the program.

For example, when you want to simulate a CPMG relaxation curve and to fit a model to the output, execute “start_cpmg_sim_fit()”.

You can change parameters such as kex, pa, J, R0, Δω, offset, CPMG period, S-spin 90 deg pulse length, I-spin 90 deg pulse length, miscalibration of S-spin RF power, as well as the kind of CPMG methods by modifying the initial part of the program (lines from 31 to 64).

The method type can be specified using the 'hpi' parameter, with the following options:

5, the 1H-CW in-phase method with pi pulse phases of {y, y, x, -x}
6, the 1H-CW in-phase method with a pi pulse at the center of the CPMG period
7, the conventional relaxation compensated CPMG method
8, the newly developed AFTAC method

The program also includes a function for fitting to relaxation-dispersion curve data that have been actually measured. Type “start_cpmg_fit(finnam='1x.inp', foutnam='1xout.txt')” with the names of the input and output files specified.

The input file (e.g. 1x.inp) looks like follows.

# comment
# cpmg_frequency_(Hz)   R2eff_(/s)   uncertainty_(any_number)
50.000000	94.018909	0.463558
100.000000	44.571446	0.065021
150.000000	25.978784	0.032330
200.000000	17.267805	0.024063
250.000000	12.571358	0.020835
300.000000	10.142765	0.019441
350.000000	8.435333	0.018562
400.000000	7.166242	0.017959
450.000000	6.555614	0.017683
50.000000	87.021661	0.350459
500.000000	5.808568	0.017357
550.000000	5.823257	0.017364
600.000000	5.347197	0.017163
650.000000	4.856996	0.016962
700.000000	4.826015	0.016950
750.000000	4.683501	0.016893
800.000000	4.536786	0.016834
