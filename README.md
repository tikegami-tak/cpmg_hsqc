# cpmg_hsqc

The programs simulate relaxation dispersion curves for NMR experiments.

The Perl and Python programs both simulate four kinds of relaxation dispersion curves. They are

1)	the conventional relaxation compensated CPMG method,
2)	the 1H-CW in-phase method with a 90deg-phase-shifted pi pulse at the center of the CPMG period,
3)	the 1H-CW in-phase method with pi pulse phases of {y, y, x, -x}
4)	the newly developed AFTAC method.

The Python program also includes a function for fitting to relaxation-dispersion curve data that have been simulated or actually measured.

The Mathematica program is the same as the python one, but has no function of fitting.

The Perl program also simulates relaxation-dispersion curves, but implements no exchange phenomenon. Instead it can deal with SI, SI2, and SI3 spin systems.

More details will be described soon.

References:

(1) T Konuma, A Nagadoi, J Kurita, & T Ikegami,
  Analysis of Artifacts Caused by Pulse Imperfections in CPMG Pulse Trains in NMR Relaxation Dispersion Experiments.
  Magnetochemistry 2018, 4(3), 33.
  https://www.mdpi.com/2312-7481/4/3/33

(2) Konuma T, Kurita J, & Ikegami T,
  Cpmg Pulse Sequence for Relaxation Dispersion that Cancels Artifacts Independently of Spin States.
  Available at SSRN:
  https://ssrn.com/abstract=4290623 or http://dx.doi.org/10.2139/ssrn.4290623
