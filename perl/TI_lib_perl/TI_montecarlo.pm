=head1 NAME

TI_montecarlo - Monte Carlo simulation

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

package TI_montecarlo;

use strict;
use Carp;

#U $x1 or ($x1, $x2) = TI_montecarlo->rand_gauss($ave, $sigma);
sub rand_gauss
{
	shift;
	my ($ave, $sigma) = @_;
	my ($pi, $r1, $r2, $fact, $x1, $x2, $ret_x1, $ret_x2);

	$pi = atan2(1, 1) * 4;

	do { $r1 = rand(); } while ($r1 <= 0);
	$r2 = rand();

	$fact = sqrt(abs(-2.0 * log($r1)));	# ln, instead of log10
	$x1 = $fact * cos(2.0 * $pi * $r2);
	$x2 = $fact * sin(2.0 * $pi * $r2);

	$ret_x2 = $sigma * $x1 + $ave;
	$ret_x1 = $sigma * $x2 + $ave;

	return wantarray ? ($ret_x1, $ret_x2) : $ret_x1;
}

#U $x1 or ($x1, $x2) = TI_montecarlo->rand_gauss2($ave, $sigma);
sub rand_gauss2
{
	shift;
	my ($ave, $sigma) = @_;
	my ($pi, $r1, $r2, $w, $fact, $x1, $x2, $ret_x1, $ret_x2);

	$pi = atan2(1, 1) * 4;

	do
	{
		$r1 = 2 * rand() - 1;
		$r2 = 2 * rand() - 1;
		$w = $r1 * $r1 + $r2 * $r2;
	} while ($w >= 1);

	$fact = sqrt((-2 * log($w)) / $w);	# ln, instead of log10
	$x2 = $fact * $r1;
	$x1 = $fact * $r2;

	$ret_x2 = $sigma * $x2 + $ave;
	$ret_x1 = $sigma * $x1 + $ave;

	return wantarray ? ($ret_x1, $ret_x2) : $ret_x1;
}

#U ($av, $sd, $num) = TI_montecarlo->ave_sigma(@data);
sub ave_sigma
{
	shift;
	my (@dat) = @_;
	my ($edat, $num, $ave, $sd);

	$ave = 0.0;
	$sd = 0.0;
	$num = 0;

	foreach $edat (@dat)
	{
		$ave += $edat;
		$num++;
	}
	if (!$num)
	{
		croak "The number of the array is 0 [ave_sigma]";
	}
	$ave /= $num;

	foreach $edat (@dat)
	{
		$sd += (($edat - $ave) * ($edat - $ave));
	}
	$sd = sqrt($sd / $num);

	return ($ave, $sd, $num);
}

#U ($eval, $numev, @par_low) = lsquare_simplex(10000, 1.0e-7, \&func, \@st_par, \@uni_par);
sub lsquare_simplex
{
	my ($self, $eval_max, $ftol, $func, $st_par, $uni_par) = @_;
	my ($ndim, $i, $j, $i_hi, $ramda, $i_nexthi, $i_low, $rtol);
	my ($swap, $ret_val, $neval, $eval_try, $eval_tmp, $cent);
	my ($parm, @evalu, @psum, @ptry);

	$ret_val = 1;
	$neval = 0;
	$parm = [[]];
	@evalu = ();
	@psum = ();
	@ptry = ();

	$ndim = scalar @{$st_par};
	if ($ndim != (scalar @{$uni_par}))
	{
		croak "The start and unit parameter numbers are not the same [lsquare_simplex]";
	}

	for ($i = 0; $i <= $ndim; $i++)
	{
		$parm->[$i] = [];
		for ($j = 0; $j < $ndim; $j++)
		{
			$ramda = rand(20) - 10;	# -10.00 ~ +10.00
			$parm->[$i][$j] = $st_par->[$j] + $ramda * $uni_par->[$j];
		}
		$evalu[$i] = $func->(@{$parm->[$i]});
	}

	for ($j = 0; $j < $ndim; $j++)
	{
		$psum[$j] = 0.0;
		for ($i = 0; $i <= $ndim; $i++) { $psum[$j] += $parm->[$i][$j]; }
	}

	for (;;)
	{
		if ($evalu[0] > $evalu[1])	{ $i_hi = 0; $i_nexthi = 1; }
		else						{ $i_hi = 1; $i_nexthi = 0; }
		$i_low = 0;

		for ($i = 0; $i <= $ndim; $i++)
		{
			if ($evalu[$i] <= $evalu[$i_low]) { $i_low = $i; }
			if ($evalu[$i] > $evalu[$i_hi]) { $i_nexthi = $i_hi; $i_hi = $i; }
			elsif (($evalu[$i] > $evalu[$i_nexthi]) && ($i != $i_hi)) { $i_nexthi = $i; }
		}

		$rtol = 2.0 * abs($evalu[$i_hi] - $evalu[$i_low]) / (abs($evalu[$i_hi]) + abs($evalu[$i_low]));
		if ($rtol < $ftol)	# END condition
		{
			$swap = $evalu[0];
			$evalu[0] = $evalu[$i_low];
			$evalu[$i_low] = $swap;

			for ($j = 0; $j < $ndim; $j++)
			{
				$swap = $parm->[0][$j];
				$parm->[0][$j] = $parm->[$i_low][$j];
				$parm->[$i_low][$j] = $swap;
			}
			$ret_val = 0;
			last;	# END
		}

		if ($neval > $eval_max)
		{
			# printf STDERR "Downhill simplex : over %d EXIT\n", $eval_max;
			$ret_val = 1;
			last;
		}

		$neval += 2;

		for ($j = 0; $j < $ndim; $j++)
		{
			$cent = ($psum[$j] - $parm->[$i_hi][$j]) / $ndim;
			$ptry[$j] = $cent -($parm->[$i_hi][$j] - $cent);
		}
		$eval_try = $func->(@ptry);

		if ($eval_try < $evalu[$i_hi])
		{
			$evalu[$i_hi] = $eval_try;
			for ($j = 0; $j < $ndim; $j++)
			{
				$psum[$j] += ($ptry[$j] - $parm->[$i_hi][$j]);
				$parm->[$i_hi][$j] = $ptry[$j];
			}
		}

		if ($eval_try <= $evalu[$i_low])
		{
			for ($j = 0; $j < $ndim; $j++)
			{
				$cent = ($psum[$j] - $parm->[$i_hi][$j]) / $ndim;
				$ptry[$j] = $cent + 2.0 * ($parm->[$i_hi][$j] - $cent);
			}
			$eval_try = $func->(@ptry);
			if ($eval_try < $evalu[$i_hi])
			{
				$evalu[$i_hi] = $eval_try;
				for ($j = 0; $j < $ndim; $j++)
				{
					$psum[$j] += ($ptry[$j] - $parm->[$i_hi][$j]);
					$parm->[$i_hi][$j] = $ptry[$j];
				}
			}
		}
		elsif ($eval_try >= $evalu[$i_nexthi])
		{
			$eval_tmp = $evalu[$i_hi];

			for ($j = 0; $j < $ndim; $j++)
			{
				$cent = ($psum[$j] - $parm->[$i_hi][$j]) / $ndim;
				$ptry[$j] = $cent + 0.5 * ($parm->[$i_hi][$j] - $cent);
			}
			$eval_try = $func->(@ptry);

			if ($eval_try < $evalu[$i_hi])
			{
				$evalu[$i_hi] = $eval_try;
				for ($j = 0; $j < $ndim; $j++)
				{
					$psum[$j] += ($ptry[$j] - $parm->[$i_hi][$j]);
					$parm->[$i_hi][$j] = $ptry[$j];
				}
			}

			if ($eval_try >= $eval_tmp)
			{
				for ($i = 0; $i <= $ndim; $i++)
				{
					if ($i != $i_low)
					{
						for ($j = 0; $j < $ndim; $j++)
						{
							$psum[$j] = ($parm->[$i][$j] + $parm->[$i_low][$j]) * 0.5;
							$parm->[$i][$j] = $psum[$j];
						}
						$evalu[$i] = $func->(@psum);
					}
				}
				$neval += $ndim;
				for ($j = 0; $j < $ndim; $j++)
				{
					$psum[$j] = 0.0;
					for ($i = 0; $i <= $ndim; $i++) { $psum[$j] += $parm->[$i][$j]; }
				}
			}
		}
		else { $neval--; }
	}

	if ($ret_val)
	{
		for ($j = 0; $j < $ndim; $j++)
		{
			$parm->[0][$j] = 0;
		}
		$evalu[0] = -1;		# must be changed to INFINITY.
		$neval = 0;
	}

	return ($evalu[0], $neval, @{$parm->[0]});
}

#U $cof_r = TI_montecarlo->linear_correlation(@data);
sub linear_correlation
{
	shift;
	my (@dat) = @_;
	my @dat_x = ();
	my @dat_y = ();
	my @word = ();
	my ($i, $edat, $num, $ave_x, $ave_y);
	my ($tmp1, $tmp2, $tmp3, $tmp4);

	$ave_x = 0.0;
	$ave_y = 0.0;
	$num = 0;

	foreach $edat (@dat)
	{
		chomp($edat);
		local $_ = $edat; @word = split;
		$dat_x[$num] = $word[0];
		$dat_y[$num] = $word[1];
		$ave_x += $word[0];
		$ave_y += $word[1];
		$num++;
	}
	if (!$num)
	{
		croak "The number of the array is 0 [linear_correlation]";
	}
	$ave_x /= $num;
	$ave_y /= $num;

	$tmp1 = 0.0;
	$tmp2 = 0.0;
	$tmp3 = 0.0;
	for ($i = 0; $i < $num; $i++)
	{
		$tmp1 += (($dat_x[$i] - $ave_x) * ($dat_y[$i] - $ave_y));
		$tmp2 += (($dat_x[$i] - $ave_x)**2);
		$tmp3 += (($dat_y[$i] - $ave_y)**2);
	}
	$tmp4 = sqrt($tmp2) * sqrt($tmp3);
	if ($tmp4 == 0)
	{
		croak "Cannot be devided by 0 [linear_correlation]";
	}

	return ($tmp1/$tmp4);
}

1;
