#! /auto/home/bs/ikegami/perl/bin/perl

$PI = atan2(1, 1) * 4;

print "\n#### Adiabatic pulse ####\n\n";

print " File name :";
$fname = <STDIN>; chomp $fname;

print " Number of splits :";
$nums = <STDIN>; chomp $nums;

print " Amplitude power index (like 40 in WURST40) :";
$powind = <STDIN>; chomp $powind;

print " Frequency sweep range (Hz) :";
$sweeprange = <STDIN>; chomp $sweeprange;								# Hz

print " Q factor on resonance :";
$q_fact = <STDIN>; chomp $q_fact;

print " Which do you choose between the pulse duration (1) or the B1 field (2) ? : ";
$ans = <STDIN>; chomp $ans;

if ($ans eq '1')
{
	print " Pulse duration (ms) :";
	$pulselen = <STDIN>; chomp $pulselen; $pulselen *= 0.001;				# sec
	$sweeprate = $sweeprange / $pulselen;									# Hz/s
	$b1field = sqrt($q_fact * $sweeprate / (2 * $PI));
}
else
{
	print " B1 field strength (Hz) :";
	$b1field = <STDIN>; chomp $b1field;										# Hz
	$sweeprate = $b1field * $b1field * 2 * $PI / $q_fact;
	$pulselen = $sweeprange / $sweeprate;
}

$durp = $pulselen / $nums;

open(SP, "> $fname") || die "Can not open the $fname file.\n";

printf SP "# WURST%d  %d points\n", $powind, $nums;
printf SP "# Pulse duration %lf ms\n", $pulselen*1000.0;
printf SP "# B1 field strength %lf kHz\n", $b1field/1000.0;
printf SP "# Frequency sweep range %lf kHz\n", $sweeprange/1000.0;
printf SP "# Q factor on resonance %.2f\n", $q_fact;
print  SP "#\n";

$amp_arg = -$PI/2;
$amp_inc = $PI / ($nums - 1);
$tim = -($nums/2.0 - 0.5) * $durp;

for ($i = 0; $i < $nums; $i++)
{
	$samp = sin($amp_arg);
	if ($samp < 0) { $samp = -$samp; }
	$amp_f = 100.0 * (1- $samp ** $powind);
	$amp_arg += $amp_inc;
	$phas_f = ($sweeprate/2.0) * $tim * $tim * 360;
	$tim += $durp;
	printf SP "%.5lf,	%.5lf\n", $amp_f, $phas_f;
#	printf SP "%d	%.5lf\n", $i, $phas_f;
}

close SP;
