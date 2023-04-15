# cpmg_hsqc

These programs simulate CPMG experiments for NMR.

Both the Perl and Python programs simulate four kinds of relaxation dispersion curves. They are conventionally used relaxation-compensated CPMG method, the 1H-CW in-phase method, the yyx(-x) in-phase CPMG method, and a newly developed AFTAC-CPMG method described in the references (1, 2).

The Python program includes a function for fitting to those relaxation curve data.

The Perl program also simulates similar relaxation-dispersion curves, but using the evolution of density matrices.

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
