# set nokey
set xrange [0:1000]
set yrange [-4:14]
set title "CPMG"
plot "mag2.dat" with lines linewidth 3, 0 with lines linewidth 1
# plot "mag.dat" with lines linewidth 3, "mag1.dat" with lines linewidth 3, 0 with lines linewidth 1
# plot "mag.dat" with errorbars
pause -1 "Hit return to continue"

plot "mag3.dat" with lines linewidth 3, 0 with lines linewidth 1
pause -1 "Hit return to continue"
