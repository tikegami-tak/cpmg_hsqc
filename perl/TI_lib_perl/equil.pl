#!/usr/sbin/perl

#-----------------------------------------------------
$init_p = 0.5e-3;		# protein conc. mol/L
$init_f = 2 * $init_p;		# ligand conc. mol/L
$Kdis   = 100e-6;		# dissociation const. mol/L
#-----------------------------------------------------

# Kdis = (init_p - x) * (init_f - x) / x

$b = $init_p + $init_f + $Kdis;
$tmp = $b**2 -4 * $init_p * $init_f;
$bac = sqrt($tmp);

$x1 = ($b + $bac)/2;
$x2 = ($b - $bac)/2;

$rat1 = $x1 / $init_p * 100;
$rat2 = $x2 / $init_p * 100;

if ($rat1 <= 100)
{
	printf "%.3f percent of the protein is in bound form.\n", $rat1;
}
if ($rat2 >= 0)
{
	printf "%.3f percent of the protein is in bound form.\n", $rat2;
}
