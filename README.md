# cpmg_hsqc

The programs simulate relaxation dispersion curves for NMR experiments. The Perl and Python programs both simulate four kinds of relaxation dispersion curves:

1. The conventional relaxation-compensated CPMG method
2. The 1H-CW in-phase method with a 90-degree phase-shifted 180-degree pulse at the center of the CPMG period
3. The 1H-CW in-phase method with 180-degree pulse phases of {y, y, x, -x}
4. The newly developed AFTAC method

The Python program also includes a function for fitting relaxation-dispersion curve data that have been simulated or actually measured.
The Mathematica program provides the same results as the Python program but does not have a fitting function.
The Perl program also simulates relaxation-dispersion curves but does not implement the exchange phenomenon. Instead, it can handle SI, SI2, and SI3 spin systems.

More details are provided in the following reference:

References:

(1) T Konuma, J Kurita, & T Ikegami,
  CPMG pulse sequence for relaxation dispersion that cancels artifacts independently of spin states.
  J. Magn. Reson. 2023, 352:107489.
  doi: 10.1016/j.jmr.2023.107489. Epub 2023 May 21.

(2) T Konuma, A Nagadoi, J Kurita, & T Ikegami,
  Analysis of Artifacts Caused by Pulse Imperfections in CPMG Pulse Trains in NMR Relaxation Dispersion Experiments.
  Magnetochemistry 2018, 4(3), 33.
  https://www.mdpi.com/2312-7481/4/3/33

(3) Konuma T, Kurita J, & Ikegami T,
  Cpmg Pulse Sequence for Relaxation Dispersion that Cancels Artifacts Independently of Spin States.
  Available at SSRN:
  https://ssrn.com/abstract=4290623 or http://dx.doi.org/10.2139/ssrn.4290623
