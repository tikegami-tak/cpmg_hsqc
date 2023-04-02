#!/usr/sbin/perl

$H_BAR = 1.05457266;			# e-34 J.s
$H_PLK = 6.626176;				# e-34 J.s
$K_CONS = 1.380662;				# e-23 J/K
$M_PI = 3.141592653;
$GYRO_1H = 26.752196e+7;
$GYRO_15N = -2.712621e+7;
$AVOGA = 6.022169;				# e+23

#-----------------------------------------------------

$nu_a = 800;					# MHz
$nu_b = $nu_a * $GYRO_15N / $GYRO_1H;
$dis_hn = 1.02;					# H-N Angstroem
$dis_hm = 30.0;					# H-M Angstroem
$temperat = 303.0;				# Temperature Kelvin
$ct = 1;						# cos(t)
$st = 0;						# sin(t)
$cdp = 1;						# cos(2f)
$xa = 0;						# Magnetic susceptibility Xa -33 m^3/mol
$xr = 0;						# Magnetic susceptibility Xr -33 m^3/mol

# The magnetic susceptibility must be sometimes devided by Avogadoro-number.
# -10 m^3/mol
# $xa = -295.0/$AVOGA;

# $xa *= (4.0 * $M_PI * 1.0e-6);	# when Xa -33 cm^3 /mol
# $xr *= (4.0 * $M_PI * 1.0e-6);	# when Xr -33 cm^3 /mol

#-----------------------------------------------------

# $xa_Hz = 0.920;				# Hz
# $xr_Hz = -0.676;				# Hz
$xa_Hz = -0.320;				# Hz
$xr_Hz = 0.209;					# Hz
$xa = $xa_Hz * 4.0*$M_PI*100 /$nu_a /$nu_b /$H_PLK *15.0 *$K_CONS *$temperat *($dis_hn**3);
$xr = $xr_Hz * 4.0*$M_PI*100 /$nu_a /$nu_b /$H_PLK *15.0 *$K_CONS *$temperat *($dis_hn**3);

$str = $xa * (3.0 * $ct * $ct - 1.0) + 1.5 * $xr * $st * $st * $cdp;

$cnst_dip = -$nu_a * $nu_b * $H_PLK / 15.0 / $K_CONS / $temperat / ($dis_hn**3);
$dip = ($str * $cnst_dip * 0.01);			# e(6+6-34+23+10+10+10-33)
$dip /= (4.0 * $M_PI);						# MKSA unit

$cnst_pcs = 1.0 / 12.0 / $M_PI / ($dis_hm**3);
$pcs = ($str * $cnst_pcs * 1000.0);			# e(-33+10+10+10+6)

print(" Xa = $xa Xr = $xr e-33 m^3/mol.\n");
print(" Residual dipolar coupling = $dip Hz.\n");
print(" Pseudo contact shift = $pcs ppm.\n");
