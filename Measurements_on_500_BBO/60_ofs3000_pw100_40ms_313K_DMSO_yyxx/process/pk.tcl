#!/bin/sh
# The next line restarts using nmrWish \
exec nmrWish "$0" -- "$@"

message .msg -relief raised -bg blue -fg white \
 -width 30c -text "Detecting Peaks ..."

pack    .msg
update

set tabName  test.tab
set specName ./Cpmg_1.ft
set tabCount 1

set tabDir [file dirname $tabName]

if {![file exists $tabDir]} {file mkdir $tabDir}


set thisSpecName $specName
set thisTabName  $tabName

set x1      171
set xN      311
set xInc    141
set xExtra  2
set xLast   [expr $xN + $xExtra + 1]

set y1      325
set yN      377
set yInc    53
set yExtra  1
set yLast   [expr $yN + $yExtra + 1]

    set yFirst  $y1

while {$yFirst <= 1 + $yN - $yExtra} \
   {
    set yNext [expr $yFirst+$yInc+2*$yExtra-1]
    if {$yNext > $yLast} {set yNext $yLast}

    set xFirst  $x1

while {$xFirst <= 1 + $xN - $xExtra} \
   {
    set xNext [expr $xFirst+$xInc+2*$xExtra-1]
    if {$xNext > $xLast} {set xNext $xLast}

    readROI -roi 1 \
       -ndim 2 -in $thisSpecName \
       -x X_AXIS $xFirst $xNext           \
       -y Y_AXIS $yFirst $yNext           \
       -verb

    pkFindROI -roi 1 \
      -sigma 16475.5 -pChi 0.0001 -plus 1.2428e+07 -minus -1.2428e+07 \
      -dx        2     1 \
      -idx       2     1 \
      -tol    7.31  2.00 \
      -hiAdj  1.80  1.20 \
      -lw     5.00 15.00 \
      -parent  -sinc -mask -out $thisTabName -verb

    set xFirst [expr 1 + $xNext - 2*$xExtra]
   }
    set yFirst [expr 1 + $yNext - 2*$yExtra]
   }

exit
