#!/usr/bin/perl

# ------------------------------------------------
# only for 13C CPMG_HSQC not for CPMG_TROSY
$Tcp = 40;		# Tcp(ms) in (Tcp/4-(1/2J)-Tcp/4)Pi(Tcp/4-(1/2J)-Tcp/4)
$L8 = 18;		# equals to NBL, L7 = 1 ... L8
$max_nu = 800;		# max = max_nu Hz

# The 1H CW powers are also calculated.
$max_cw = 13000;		# The max power (Hz) for the 1H CW
# ------------------------------------------------

$maxL6 = int($max_nu * $Tcp /2000 + 0.5);
$min_nu = 2000/$Tcp;
$step_nu = ($max_nu-$min_nu)/($L8-2);

$numPi = 4*$maxL6 + 2;
print ";\n";

for ($i = 1; $i <= $L8; $i++)
{
	printf "if \"l7==%d\" goto v%02d", $i, $i;
	if ($i == 1) { print "		; reference"; }
	print "\n";

	# Equal step
	$L6array[$i] = int(($min_nu + $step_nu * ($i-2)) * $Tcp/2000 +0.5);
}
$L6array[1] = 0;

print ";\n; otherwise goto v01\n;\n";

for ($i = 1; $i <= $L8; $i++)
{
	printf "v%02d, d11\n", $i;

	printf "   \"l6=%d\"", $L6array[$i];
	printf "		; %.2f Hz", $L6array[$i]*2000/$Tcp;
	if ($i == 2) { print " v=L6/[Tcp/2]"; }
	print "\n";

	if ($i == 1)
	{
		printf "   \"d17=1m\"		; dummy";
	}
	else
	{
		$tmp = int(100000*$Tcp/8/$L6array[$i] +0.5);
		$tmp /= 100000.0;
		if ($tmp < 1.0)
		{
			printf "   \"d17=%.2fu-p21\"", $tmp*1000;
		}
		else
		{
			printf "   \"d17=%.5fm-p21\"", $tmp;
		}
	}
	if ($i == 2) { print "	; Tcp/(8*L6)"; }
	print "\n";

	$L5 = $numPi-4*$L6array[$i];
	printf "   \"l5=%d\"", $L5;
	if ($i == 2) { printf "		; L5=%d-4*L6", $numPi; }
	print "\n";

	printf "   \"d19=%.3fm-p21\"", 1000/2/$L5;
	if ($i == 2) { print "	; 1s/(2*L5)"; }
	print "\n";

	print "   d11\n";
	print "   goto v70\n";
}
print "\n;---------\n";
printf "; d14 = %d m\n", $Tcp;
printf "; l8 = %d\n", $L8;

# 1H cw powers

print "\n; ----- The following lines should be put before ze. -----\n\n";

$plw = 12;
for ($i = 2; $i <= $L8; $i++)
{
	my $k = 1;
	$p90cw = 1000.0*$Tcp/16/$L6array[$i];
	while (($p90cw/$k) > (1e+6)/$max_cw/4.0) { $k++; }
	$p90cw /= $k;
	printf "\"plw%d=plw11*pow(p11/%.2f, 2)\"", $plw, $p90cw;
	print "\n";
	$plw++;
}

print "\n; ----- The following lines should be put before every CPMG pulse train. -----\n\n";

for ($i = 1; $i <= $L8; $i++)
{
	printf "if \"l7==%d\" goto w%02d\n", $i, $i;
}

$plw = 11;
for ($i = 1; $i <= $L8; $i++)
{
	printf "w%02d, 20u pl%d:f1\n", $i, $plw;
	print "   goto w70\n";
	$plw++;
}
