#!/usr/bin/perl

# Takahisa Ikegami
# Draw the point at x=50,y=100 in the file (378*256).
# When there is no assignment

$xt = 512;
$yt = 513;
$x_box = 6;	# x_box *2 +1
$y_box = 1;	# y_box *2 +1

$st = 0;
$k = 0;

open(SHIFT, "< ../process/230424_test.tab")	|| die "Can not open the test.tab file.\n";
open(AMNO, "> extpoint.dat")	|| die "Can not open the output file for write.\n";

# 3   187.674   103.330  0.073  0.004  127.520   10.504  5170.924  4203.001   9.070

LINET:while (<SHIFT>)
{
	$line = $_;
	next LINET if (/^#/ || /^!/ || /^\n/);
	chomp;
	@word = split;
	
	if ($word[0] eq 'FORMAT')
	{
		$st = 1;
		goto LINET;
	}
	if (!$st) { goto LINET; }

	if ($word[21] ne '1')
	{
		goto LINET;
	}

	$bestX = $word[5];
	$bestY = $word[6];

	$resnum = $word[0];
	$resnam = 'X';
	$resnam3 = 'Xxx';

	$x = int($word[1] + 0.5);
	$y = int($word[2] + 0.5);

	if (($x < $x_box+1) || ($x > $xt-$x_box-1) || ($y < $y_box+1) || ($y > $yt-$y_box-1))
	{
		die "No. $resnam $resnum is out of the range. [$x, $y pnt] [$bestX, $bestY ppm].\n";
	}
	$pti = $xt * ($y-1) + $x;
	$tot_pt[$k++] = $pti;

	printf AMNO "RES	%d	%.6f	%.6f	%d	%d	%d	%s\n",
		$resnum, $bestX, $bestY, $x, $y, $pti, $resnam;

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
}
close SHIFT;
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
		die " Cannot read the input data ## error.\n";
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
