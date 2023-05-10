#!/bin/csh

set pdir = "../59_ofs3000_pw90_40ms_313K_DMSO_yyxx"

mkdir process pick_area

cp $pdir/sep*pl ./

cp $pdir/process/*com process/
cp $pdir/process/*pl process/

cp $pdir/pick_area/*pl pick_area/
cp $pdir/pick_area/indoxylsulfate_shift.txt pick_area/
cp $pdir/pick_area/indoxylsulfate_shift_fold2.txt pick_area/
