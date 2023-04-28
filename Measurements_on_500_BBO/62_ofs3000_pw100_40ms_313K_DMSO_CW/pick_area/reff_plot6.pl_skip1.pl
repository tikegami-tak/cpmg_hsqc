#!/usr/bin/perl

$seed = 463;
$title = 'DMTCA CPMG_methyl (500MHz_new_version) atom = ';
$x_box = 6;
$y_box = 4;
$box_size = ($x_box *2 +1) * ($y_box *2 +1);
$uncert0 = 14830*sqrt($box_size);	# reference
$uncert1 = 14830*sqrt($box_size);
$const_time = 0.08;	# 40 ms

# $sername[0] = "cpmg_800_xxxx" . "$i" . ".ft";

$sername[0] = "Cpmg_1.ft";
# $sername[1] = "Cpmg_2.ft";
$sername[1] = "Cpmg_3.ft";
$sername[2] = "Cpmg_4.ft";
$sername[3] = "Cpmg_5.ft";
$sername[4] = "Cpmg_6.ft";
$sername[5] = "Cpmg_7.ft";
$sername[6] = "Cpmg_8.ft";
$sername[7] = "Cpmg_9.ft";
$sername[8] = "Cpmg_10.ft";
$sername[9] = "Cpmg_11.ft";
$sername[10] = "Cpmg_12.ft";
$sername[11] = "Cpmg_13.ft";
$sername[12] = "Cpmg_14.ft";
$sername[13] = "Cpmg_15.ft";
$sername[14] = "Cpmg_16.ft";
$sername[15] = "Cpmg_17.ft";
$sername[16] = "Cpmg_18.ft";
$sername[17] = "Cpmg_19.ft";
$sername[18] = "Cpmg_20.ft";
$sername[19] = "Cpmg_21.ft";
$sername[20] = "Cpmg_22.ft";
$sername[21] = "Cpmg_23.ft";
$sername[22] = "Cpmg_24.ft";
$sername[23] = "Cpmg_25.ft";
$sername[24] = "Cpmg_26.ft";
$sername[25] = "Cpmg_27.ft";
$sername[26] = "Cpmg_28.ft";
$sername[27] = "Cpmg_29.ft";
$sername[28] = "Cpmg_30.ft";
$sername[29] = "Cpmg_31.ft";
$sername[30] = "Cpmg_32.ft";
$sername[31] = "Cpmg_33.ft";

$tim[0] = 0;            # Hz
# $tim[1] = 25;
$tim[1] = 50;
$tim[2] = 75;
$tim[3] = 100;
$tim[4] = 125;
$tim[5] = 150;
$tim[6] = 175;
$tim[7] = 200;
$tim[8] = 225;
$tim[9] = 250;
$tim[10] = 275;
$tim[11] = 300;
$tim[12] = 325;
$tim[13] = 350;
$tim[14] = 375;
$tim[15] = 400;
$tim[16] = 425;
$tim[17] = 450;
$tim[18] = 475;
$tim[19] = 500;
$tim[20] = 525;
$tim[21] = 550;
$tim[22] = 575;
$tim[23] = 600;
$tim[24] = 625;
$tim[25] = 650;
$tim[26] = 675;
$tim[27] = 700;
$tim[28] = 725;
$tim[29] = 750;
$tim[30] = 775;
$tim[31] = 800;

$ntim = 32;

if (!(-d "./reff_plot"))
{
	mkdir("reff_plot", 0777) || die "Can not make the directory reff_plot.\n";
}
system("rm reff_plot/*");

open(EXT, "< ./extpoint.dat") || die "Can not open the file extpoint.dat.\n";
$i = -1;
$sem = 0;
$k = $box_size;
LINE:while(<EXT>)
{
	next LINE if (/^#/ || /^ / || /^\n/);
	chomp;
	@word = split;
	if ($word[0] eq 'RES')
	{
		if ($k != $box_size)
		{
			die "$k != box_size ($box_size) in amino acid $i\n";
		}
		$i++;
		$lab[$i] = $word[7];
		$res[$i] = $word[1];
		$alab[$i] = "$word[1]" . "$word[7]";
		$chemx[$i] = $word[3];
		$chemy[$i] = $word[2];
		$sem = 1;
		$k = 0;
		goto LINE;
	}
	else
	{
		if ($sem == 1)
		{
			$pt[$i][$k++] = $word[0];
		}
	}
}
close(EXT);
$anum = $i + 1;

if ($ntim != scalar(@sername))
{
	die "The numbers of the times and the sers are inconsistent.\n";
}

$t = 0;
foreach $sern (@sername)
{
	open(INP, "< ../process/$sern") || die "Can not open the DAT file of $sern for read.\n";

	for ($j = 0; $j < $anum; $j++)
	{
		$alabpck[$j][$t] = 0;

		for ($k = 0; $k < $box_size; $k++)
		{
			$nbyte = 2048 + 4*($pt[$j][$k] - 1);

			seek INP, $nbyte, 0;
			$rebyte = sysread(INP, $iidata, 4);
			if ((!defined $rebyte) || ($rebyte != 4))
			{
				die " Can not read the input data of $sern. ## error.\n";
			}
			$alabpck[$j][$t] += (unpack "f", $iidata);
		}
	}
	close INP;
	$t++;
}

# print "anum:  $anum\n";

for ($j = 0; $j < $anum; $j++)
{
	$maxyplot = 0;

	$inpfile = "$alab[$j]".'.inp';
	open(INP, "> reff_plot/$inpfile") || die "Can not open the file $inpfile.\n";

	print INP "\# $title $alab[$j]\n";
	# $nseed = int(rand($seed));
	# print INP "$ntim\t$nseed\n";
	# print INP "data\n";

	# if ($alabpck[$j][0] > 80000000) { $maxyplot = 1; }

	for ($t = 1; $t < $ntim; $t++)
	{
		$reff = -log($alabpck[$j][$t] / $alabpck[$j][0])/$const_time;

		$uncert = sqrt((($uncert1/$alabpck[$j][$t])**2) + (($uncert0/$alabpck[$j][0])**2))/$const_time;

		printf INP "%f\t%f\t%f\n", $tim[$t], $reff, $uncert;
	}
	close(INP);

	$plotoutfile = "$alab[$j]" . '.plot';
	open(PLOTOUT, "> reff_plot/$plotoutfile") || die "Can not open the file $plotoutfile.\n";


	print PLOTOUT ("set nokey\n");
	print PLOTOUT ("set xrange [-10:810]\n");

	# if (!$maxyplot) { print PLOTOUT ("set yrange [-1:8]\n"); }

	print PLOTOUT "set title \"$title $alab[$j] $chemx[$j] $chemy[$j]\"\n";
	print PLOTOUT ("plot \"$inpfile\" with errorbars, \"$inpfile\" with impulses, 0 with lines\n");
	print PLOTOUT ("pause -1 \"Hit return to continue\"\n");

	close(PLOTOUT);
}
