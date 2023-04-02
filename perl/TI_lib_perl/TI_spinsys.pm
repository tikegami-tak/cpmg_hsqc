=head1 NAME

TI_spinsys - Spin system calculations

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

use strict;
use Carp;
use Math::Complex;

my $THRESH_Z = 1.0e-7;	# You can change this value.

#U $sigma = $spin_sys->make_density_matrix(3*(4-2*i), "e e", 2, "x y");
sub make_density_matrix
{
	my $self = shift;
	my @state_array = @_;
	my ($state_ele, $state, $mat, $elem, $cons);
	my ($state_mat, $retmat, $spinn);
	my ($re, $im);
	my @state_ele;

	$retmat = undef;
	$re = 1;	$im = undef;
	foreach $state (@state_array)
	{
		if ($state =~ /\d/)
		{
			if (($re != 1) || (defined $im))
			{
				croak "The number appeared twice sequentially [make_density_matrix]";
			}
			if ($state =~ /i/)
			{
				$re = ($state->Re); $im = ($state->Im);
			}
			else
			{
				$re = $state;
			}
		}
		else
		{
			#		$state =~ s/^\s+//;	$state =~ s/\s+$//;
			#		@state_ele = split(/\s+/, $state);
			local $_ = $state; @state_ele = split;

			$state_mat = undef;
			$spinn = 0;

			foreach $elem (@state_ele)
			{
				if ($elem eq 'e')
				{
					$mat = TI_matrix->new([1, 0], [0, 1]);	# E
				}
				elsif ($elem eq 'x')
				{
					$mat = TI_matrix->new([0, 0.5], [0.5, 0]);	# X
				}
				elsif ($elem eq 'y')
				{
					$mat = TI_matrix->new([0, "0 -0.5"], ["0 0.5", 0]);	# Y
				}
				elsif ($elem eq 'z')
				{
					$mat = TI_matrix->new([0.5, 0], [0, -0.5]);	# Z
				}
				elsif ($elem eq 'p')
				{
					# $cons = -1.0/sqrt(2.0);
					$cons = 1.0;
					$mat = TI_matrix->new([0, $cons], [0, 0]);	# +
				}
				elsif ($elem eq 'n')
				{
					# $cons = 1.0/sqrt(2.0);
					$cons = 1.0;
					$mat = TI_matrix->new([0, 0], [$cons, 0]);	# -
				}
				elsif ($elem eq 'a') # 1/2+Z
				{
					$mat = TI_matrix->new([1, 0], [0, 0]);	# alpha
				}
				elsif ($elem eq 'b') # 1/2-Z
				{
					$mat = TI_matrix->new([0, 0], [0, 1]);	# beta
				}
				else
				{
					croak "The symbol $elem is not defined [make_density_matrix]";
				}
	
				$spinn++;
				if (defined $state_mat)
				{
					$state_mat = $state_mat->direct_product($mat);
				}
				else
				{
					$state_mat = $mat;
				}
			}

			if (scalar(@{$self->{spin_array}}) != $spinn)
			{
				croak "The number of nuclei is inconsistent with spin system [make_density_matrix]";
			}

			if (defined $im)
			{
				$state_mat = $state_mat->scalar_multiply("$re $im");
			}
			elsif ($re != 1)
			{
				$state_mat = $state_mat->scalar_multiply($re);
			}

			if (defined $retmat)
			{
				$retmat = $retmat->add($state_mat);
			}
			else
			{
				$retmat = $state_mat;
			}

			$re = 1;	$im = undef;
		}
	}

	unless (defined $retmat)
	{
		croak "There is not a valuable argument [make_density_matrix]";
	}

	return $retmat;
}

#U $res = $spin_sys->evolve_density_matrix($sigma, $hamilt, $t);
sub evolve_density_matrix
{
	my ($self, $sigma, $hamil, $dt) = @_;
	my @res = ();

	@res = $self->sequential_evolution($sigma, $hamil, $dt, 1);
	return $res[1];
}

#U @res = $spin_sys->sequential_evolution($sigma, $hamilt, $dt, $n);
# @res has N+1 data ranging from 0 to N.
# Hamiltonian should be time-independent
# and does not contain relaxation or time-dependent RF pulses.
sub sequential_evolution
{
	my $func = (caller(0))[3];
	if (@_ != 5)
	{
		print STDERR "USAGE: \$result = \$spin_system->$func(\$sigma, \$hamilt, \$dt, \$n);\n\n";
		croak "Wrong usage";
	}

	my ($self, $sigma, $hamil, $dt, $n) = @_;
	my ($i, $tmp, $dexp, $uexp, $mdt, $r_uexp);
	my ($nrow, $super_ham, $col_dmat);
	my @res = ();

	$mdt = -"$dt";
	$res[0] = $sigma->copy_matrix;

	if ($hamil->is_hermitian_matrix)
	{
		# Since H is hermitian, exp(-i H t) is unitary matrix.
		# Hermitian-conjugate of unitary is the same as its inverse.
		$uexp = $hamil->make_identity_matrix;
		$tmp = $hamil->scalar_multiply("0 $mdt");
		$dexp = $tmp->exponential_operator;
		foreach $i (1..$n)
		{
			$uexp = $uexp->multiply($dexp);
			$r_uexp = $uexp->hermitian_conjugate;
			$tmp = $uexp->multiply($sigma);
			$res[$i] = $tmp->multiply($r_uexp);
		}
	}
	else
	{
		print STDERR "Hamiltonian is not a hermitian matrix.\n";
		$nrow = 2**(scalar(@{$self->{spin_array}}));
		$super_ham = $hamil->commutator_superoperator;
		$uexp = $super_ham->make_identity_matrix;
		$tmp = $super_ham->scalar_multiply("0 $mdt");
		$dexp = $tmp->exponential_operator;
		$col_dmat = $sigma->mat_to_column_vector;
		foreach $i (1..$n)
		{
			$uexp = $uexp->multiply($dexp);
			$tmp = $uexp->multiply($col_dmat);
			$res[$i] = $tmp->column_vector_to_mat($nrow, $nrow);
		}
	}
	return @res;
}

#U $data = $spin_sys->sequential_detect($a_sigma, $obs, $n);
# $a_sigma is a reference of an array of density matrices.
sub sequential_detect
{
	my $func = (caller(0))[3];
	if (@_ != 4)
	{
		croak "Wrong usage in $func";
	}

	my ($self, $a_sigma, $obs, $n) = @_;
	my ($tmp, $i, $re, $im);
	my @res = ();

	foreach $i (0..$n)
	{
		$tmp = ($a_sigma->[$i])->multiply($obs);
		($re, $im) = $tmp->trace;
		$res[$i] = "$re $im";
	}
	return @res;
}

#U $data = $spin_sys->detect_save_ff($sigma, $hamil, $obs, $dt, $n, $file, $rg, $relax);
sub detect_save_ff
{
	my $self = shift;
	my @ares = ();
	my @ri_data = ();
	my ($i, $rel, $re, $im, $ffdata);

	if (@_ != 8)
	{
		print STDERR "Arguments should be \$sigma, \$hamil, \$obs, \$dt, \$n, \$file, \$rg, \$relax\n";
		croak "Wrong usage";
	}
	my ($sigma, $hamil, $obs, $dt, $n, $file, $rg, $relax) = @_;

	@ares = $self->sequential_evolution($sigma, $hamil, $dt, $n);
	@ri_data = $self->sequential_detect(\@ares, $obs, $n);

	open(SER, ">> $file")   || croak "Can not open the output file $file.\n";

	for ($i = 0; $i < $n; $i++)
	{
		if (abs($relax) < $THRESH_Z)
		{
			$rel = 0;
		}
		else
		{
			$rel = (-1) * $dt * $i / $relax;
		}
		local $_ = $ri_data[$i];
		($re, $im) = split;
		$re *= ($rg * exp($rel));
		$im *= ($rg * exp($rel));
		$ffdata = pack "ff", $re, $im;
		syswrite(SER, $ffdata, 8);
	}

	close SER;
}

#U $spin_sys->product_operator_deconv($sigma);
sub product_operator_deconv
{
	my ($self, $sigma) = @_;
	my ($obs, $spinn, $obs_spin, $i, $j);
	my ($reu, $imu, $re_norm, $im_norm);
	my @ind = ();

	my @cartesian = ('e', 'x', 'y', 'z');

	$spinn = scalar(@{$self->{spin_array}});
	if ($spinn < 1)
	{
		croak "The number of the spins is less than 1 [product_operator_deconv]";
	}


	for ($i = 0; $i <= $spinn; $i++) { $ind[$i] = 0; }

	while ($ind[$spinn] < 1)
	{
		$obs_spin = '';
		for ($i = 0; $i < $spinn; $i++)
		{
			$j = $ind[$i];
			$obs_spin = "$obs_spin" . "$cartesian[$j] ";
		}
		chop $obs_spin;

		$obs = $self->make_density_matrix($obs_spin);
		($reu, $imu, $re_norm, $im_norm) = $self->vector_projection($obs, $sigma);

		if ((abs($re_norm) > $THRESH_Z) || (abs($im_norm) > $THRESH_Z))
		{
			printf("  %s   %.3f %.3f I\n", $obs_spin, $re_norm, $im_norm);
		}

		$ind[0]++;
		for ($i = 0; $i < $spinn; $i++)
		{
			if ($ind[$i] == @cartesian)
			{
				$ind[$i] = 0;
				$ind["$i"+1]++;
			}
		}
	}
}

#U $res = $spin_sys->rotate_density_matrix($sigma, "z e", $PI/2);
sub rotate_density_matrix
{
	my ($self, $sigma, $dens, $ang) = @_;
	my ($res, $dens_m);

	$dens_m = $self->make_density_matrix($dens);
	$res = $self->evolve_density_matrix($sigma, $dens_m, $ang);
	return $res;
}

1;
