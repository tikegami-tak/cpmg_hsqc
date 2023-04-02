#!/usr/bin/perl -I../../TI_lib_perl

#!\C:\Strawberry\perl\bin\perl.EXE -I\Users\tiik\Desktop\GoogleDrive\main_data_2\for_desktop\program\tisim\TI_lib_perl

# use diagnostics;
use TI_SIM;

$| = 1;		# STDOUT flush

$PI = atan2(1, 1) * 4;

# Parameters

$ph = 0;	# 0: correct, 1: wrong
$sem = 3;	# 1, n=0,1,40 Rex vs offset for the new sequence, AFTAC
			# 4, n=0,1,40 for the normal RC (even numbered pulses)
			# 2, dispersion curve Rex vs Vcpmg for the new sequence, AFTAC
			# 3, dispersion curve Rex vs Vcpmg for the normal RC (even numbered pulses)
$pmis = 0.9;
$dmis = 1.0;
$c90 = 20e-6;
$i90 = 1e-12;	# infinite power actually
$Cpwr = 2.0*$PI/($c90*4) *$pmis;
$Ipwr = 2.0*$PI/($i90*4);

$j_val = 140;
$Jcoup = $j_val*2*$PI;

$tex = 40e-3;
$kmax = 21;
$cpn = 1;
$d17 = $tex/$cpn/8.0 - $c90;
$dc = 2.2;		# just determined by simulation

$totn = 300;                          # number of points
$maxofs = 3000;
$offs = 4000;

open(MF, "> mag3.dat") || die "Can not open the output file.\n";

$a = TI_SIM->new("C K L M");	# K, L, and M are equivalent 1H protons in CH3.

# $obs3 = $obs1->add($obs2);
# $obs3 = $obs1->subtract($obs2);
# $obs = $obs3->scalar_multiply(0.5);

if (($sem == 3) || ($sem == 4))
{
	$obs = $a->make_density_matrix(
		4/sqrt(3), "x e z z",
		4/sqrt(3), "x z e z",
		4/sqrt(3), "x z z e");
	# $obs = $a->make_density_matrix(4/sqrt(3), "x a z z", 4/sqrt(3), "x z a z", 4/sqrt(3), "x z z a");
	# $obs = $a->make_density_matrix(4/sqrt(3), "x b z z", 4/sqrt(3), "x z b z", 4/sqrt(3), "x z z b");
	# $obs = $a->make_density_matrix(4/sqrt(3), "x a z z", 4/sqrt(3), "x z a z", 4/sqrt(3), "x z z a");
}
else
{
	$obs = $a->make_density_matrix(
		2/sqrt(3), "y z e e",
		2/sqrt(3), "y e z e",
		2/sqrt(3), "y e e z");
}

if ($sem == 1)
{
	for ($j = 0; $j < $totn; $j++)
	{
		printf "%d ", $j-$totn/2.0;

		$offs = ($j-$totn/2.0)*2*$maxofs/$totn;

		$hamil = $a->make_density_matrix(
			2*$PI*$offs, "z e e e",
			$Jcoup, "z z e e",
			$Jcoup, "z e z e",
			$Jcoup, "z e e z");

		$renx = &cpmg_pulse_aftac;
		printf MF "%.3f %.6f\n", $offs, $renx;
	}
}
elsif ($sem == 4)
{
	for ($j = 0; $j < $totn; $j++)
	{
		printf "%d ", $j-$totn/2.0;

		$offs = ($j-$totn/2.0)*2*$maxofs/$totn;

		$hamil = $a->make_density_matrix(
			2*$PI*$offs, "z e e e",
			$Jcoup, "z z e e",
			$Jcoup, "z e z e",
			$Jcoup, "z e e z");

		$renx = &cpmg_pulse_rc_conventional;
		printf MF "%.3f %.6f\n", $offs, $renx;
	}
}
elsif ($sem == 2)
{
	$hamil = $a->make_density_matrix(
		2*$PI*$offs, "z e e e",
		$Jcoup, "z z e e",
		$Jcoup, "z e z e",
		$Jcoup, "z e e z");

	for ($cpn = 0; $cpn < $kmax; $cpn++)
	{
		print "$cpn ";

		if ($cpn == 0)
		{
			$ren0 = &cpmg_pulse_aftac;
		}
		else
		{
			$d17 = $tex/$cpn/8.0 - $c90;
			$cpv = 2.0*$cpn/$tex;
			$renx = &cpmg_pulse_aftac;
			$rex = -(log($renx/$ren0))/$tex;
			printf MF "%.3f %.6f\n", $cpv, $rex;
		}
	}
}
elsif ($sem == 3)
{
	$hamil = $a->make_density_matrix(
		2*$PI*$offs, "z e e e",
		$Jcoup, "z z e e",
		$Jcoup, "z e z e",
		$Jcoup, "z e e z");

	for ($cpn = 0; $cpn < $kmax; $cpn++)
	{
		print "$cpn ";

		if ($cpn == 0)
		{
			$ren0 = &cpmg_pulse_rc_conventional;
		}
		else
		{
			$d17 = $tex/$cpn/8.0 - $c90;
			$cpv = 2.0*$cpn/$tex;
			$renx = &cpmg_pulse_rc_conventional;
			$rex = -(log($renx/$ren0))/$tex;
			printf MF "%.3f %.6f\n", $cpv, $rex;
		}
	}
}
close MF;
print "\n";

sub cpmg_pulse_aftac
{
	my $k, $i;

	for ($k = 0; $k < 2; $k++)
	{
		$sigma = $a->make_density_matrix(
			2/sqrt(3), "y z e e",
			2/sqrt(3), "y e z e",
			2/sqrt(3), "y e e z");

		for ($i = 0; $i < $cpn; $i++)	# the 1st CPMG
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*2);	# 1 y
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		# 0 x
		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*60/90);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*150/90);
		$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*150/90);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*60/90);
		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));

		for ($i = 0; $i < $cpn; $i++)	# the 2nd CPMG
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*2);	# 0 x
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		if ($ph == 0)	# 0/2 x/-x correct
		{
			if ($k == 0)
			{
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*300/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI, $c90*60/90);
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
			}
			else
			{
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI, $c90*300/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
			}
		}
		else	# 1/3 y/-y wrong
		{
			if ($k == 0)
			{
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*300/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*60/90);
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
			}
			else
			{
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*300/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*60/90);
				$sigma = $a->evolve_density_matrix($sigma, $hamil, 1.22/1000);
			}
		}

		for ($i = 0; $i < $cpn; $i++)	# the 3rd CPMG
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*2);	# 2 -x
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);

		}

		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*150/90);
		$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*150/90);
		$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));

		for ($i = 0; $i < $cpn; $i++)	# the 4th CPMG
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*2);	# 1 y
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		($reu, $imu, $re_norm, $im_norm) = $a->vector_projection($obs, $sigma);

		if ($k == 0)
		{
			$ren = $re_norm;
		}
		else
		{
			$ren += $re_norm;
			$ren *= 0.5;
		}
	}
	return $ren;
}

sub cpmg_pulse_rc_conventional
{
	my $k, $i;

	for ($k = 0; $k < 2; $k++)
	{
		$sigma = $a->make_density_matrix(
			2/sqrt(3), "y z e e",
			2/sqrt(3), "y e z e",
			2/sqrt(3), "y e e z");

		for ($i = 0; $i < $cpn; $i++)
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*2);
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		for ($i = 0; $i < $cpn; $i++)
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*2);
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));

		if ($ph == 0)	# right
		{
			if ($k == 0)
			{
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*60/90);
			}
			else
			{
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*60/90);
			}
		}
		else	# wrong
		{
			if ($k == 0)
			{
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*60/90);
			}
			else
			{
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*60/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "K", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "L", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "M", $hamil, $Ipwr, 0, $i90*2);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, -$PI/2, $c90*150/90);
				$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, $PI/2, $c90*60/90);
			}
		}
		$sigma = $a->evolve_density_matrix($sigma, $hamil, $dmis*(1/4/$j_val-$c90*2.0/$PI*$dc));

		for ($i = 0; $i < $cpn; $i++)
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*2);
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		for ($i = 0; $i < $cpn; $i++)
		{
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
			$sigma = $a->apply_pulse($sigma, "C", $hamil, $Cpwr, 0, $c90*2);
			$sigma = $a->evolve_density_matrix($sigma, $hamil, $d17);
		}

		($reu, $imu, $re_norm, $im_norm) = $a->vector_projection($obs, $sigma);

		if ($k == 0)
		{
			$ren = $re_norm;
		}
		else
		{
			$ren += $re_norm;
			$ren *= 0.5;
		}
	}
	# print "R $ren   ";
	return $ren;
}
