#!/usr/bin/perl

# Takahisa Ikegami
# Draw the point at x=50,y=100 in the file (378*256).

$xt = 512;
$yt = 206;
$x_pt = 2512.562988 / (500.13*0.251449530) / $xt;	# 0.101329118
$y_pt = 1005.859375 / 500.13 / $yt;
$x_ct = 120.179;	# 122.679;
$y_ct = (6.5+8.5)/2.0;
$limX = 0.2;
$limY = 0.03;
$devchemX = 0; # 93/2/500.13/0.101329118;	# 1J(HN)/2 for 15N
$devchemY = 0; # -93/2/500.13;		# 1J(HN)/2 for 1H
$x_box = 5;	# x_box *2 +1
$y_box = 2;	# y_box *2 +1

$a = 0;
$st = 0;
open(SHIFT, "< ../process/190411_test.tab")	|| die "Can not open the test.tab file.\n";

# 3   187.674   103.330  0.073  0.004  127.520   10.504  5170.924  4203.001   9.070

LINET:while (<SHIFT>)
{
	$line = $_;
	next LINET if (/^#/ || /^!/ || /^\n/ || /^\r/);
	chomp;
	@word = split;
	
	if ($word[0] eq 'FORMAT')
	{
		$st = 1;
		goto LINET;
	}
	if (!$st) { goto LINET; }

	$chemX[$a] = $word[5];
	$chemY[$a] = $word[6];
	$numX[$a] = $word[1];
	$numY[$a] = $word[2];
	$a++;
}
close SHIFT;

# open(AMN, "< indoxylsulfate_shift_fold2.txt") || die "Can not open the chemical shift file for read.\n";
open(AMN, "< indoxylsulfate_shift.txt") || die "Can not open the chemical shift file for read.\n";
open(AMNO, "> extpoint.dat")	|| die "Can not open the output file for write.\n";

# e3	120.108    7.703
# 33   Ala     A	127.456	9.239

$i = 0;
$k = 0;
LINE:while (<AMN>)
{
	$line = $_;
	next LINE if (/^#/ || /^!/ || /^\n/ || /^\r/);
	chomp;
	@word = split;
	
#	$_ = $word[0];
#	/^(\w)(\d+)N-H$/;
#	$resnum = $2;
#	$resnam = $1;
#	/^<(\w)(\d+).N,HN>$/;

	$resnum = $word[0];
	$resnam = $word[0];
	$chm_x = $word[1] + $devchemX;
	$chm_y = $word[2] + $devchemY;
	$resnam3 = $word[0];

	$minvec = 1000;
	for ($j = 0; $j < $a; $j++)
	{
		$difX = abs($chemX[$j] - $chm_x);
		$difY = abs($chemY[$j] - $chm_y);
		if (($difX < $limX) && ($difY < $limY))
		{
			$vec = ($difX / $limX)**2 + ($difY / $limY)**2;
			$vec = sqrt($vec);
			if ($vec < $minvec)
			{
				$minvec = $vec;
				$x = int($numX[$j] + 0.5);
				$y = int($numY[$j] + 0.5);
				$bestX = $chemX[$j];
				$bestY = $chemY[$j];
			}
		}
	}

	if ($minvec == 1000)
	{
		$x = int(($x_ct - $chm_x) / $x_pt + $xt / 2 + 0.5);
		$x += 1;	# Adjustment
		$y = int(($y_ct - $chm_y) / $y_pt + $yt / 2 + 0.5);
		$y += 1;	# Adjustment
		$bestX = $chm_x;
		$bestY = $chm_y;
		printf "%s %d %.6f %.6f is not found in test.tab\n", $resnam, $resnum, $chm_x, $chm_y;
	}

	if (($x < $x_box+1) || ($x > $xt-$x_box-1) || ($y < $y_box+1) || ($y > $yt-$y_box-1))
	{
		die "No. $resnam $resnum is out of the range. [$x, $y].\n";
	}
	$pt[$i] = $xt * ($y-1) + $x;
	$tot_pt[$k++] = $pt[$i];

	printf AMNO "RES	%d	%.6f	%.6f	%d	%d	%d	%s\n",
		$resnum, $bestX, $bestY, $x, $y, $pt[$i], $resnam;

	for ($xp = -$x_box; $xp <= $x_box; $xp++)
	{
		for ($yp = -$y_box; $yp <= $y_box; $yp++)
		{
			$pt1 = $xt * ($y -1 + $yp) + $x + $xp;
			printf AMNO "	%d\n", $pt1;

			if (($xp == -$x_box) || ($xp == $x_box) || ($yp == -$y_box) || ($yp == $y_box))
			{
				$tot_pt[$k++] = $pt1;
			}
		}
	}

	$i++;
}
close AMN;
close AMNO;

open(INP, "< ../process/Cpmg_1.ft")	|| die "Can not open the input DAT file for read.\n";
open(OUTP, "> temp.ft")	|| die "Can not open the output file for write.\n";

$rebyte = sysread(INP, $iidata, 2048);
if ((!defined $rebyte) || ($rebyte != 2048))
{
	die " Cannot read the header data ## error.\n";
}
$wrbyte = syswrite(OUTP, $iidata, 2048);

$cnt = 1;
while (1)
{
	$rebyte = sysread(INP, $iidata, 4);
	if (!defined $rebyte)
	{
		die " Can not read the input data ## error.\n";
	}
	elsif ($rebyte == 0)
	{
		last;
	}
	$val = unpack "f", $iidata;

	foreach $ept (@tot_pt)
	{
		if ($cnt == $ept)
		{
			if ($val > 0)
			{
				$val = -50000000;
			}
			else
			{
				$val = 50000000;
			}
			$iidata = pack "f", $val;
		}
	}

	$wrbyte = syswrite(OUTP, $iidata, 4);
	$cnt++;
}

close INP; close OUTP;
