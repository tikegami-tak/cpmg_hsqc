set nokey
set xrange [-50:2100]
set title "DMSO CPMG_methyl (500MHz_cw313K_version) atom =  1X 2.545000 43.697000"
plot "1X.inp" with errorbars, "1X.inp" with impulses, 0 with lines
pause -1 "Hit return to continue"
