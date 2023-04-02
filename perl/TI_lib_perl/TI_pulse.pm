=head1 NAME

TI_pulse - Pulse calculation

=head1 AUTHOR

Takahisa Ikegami <ikegamit@yokohama-cu.ac.jp>

=cut

use strict;
use Carp;

#U $res = $spin_sys->apply_pulse($sigma, "I", $hamil, 2*$PI*3000, $PI/2, 10.0e-6);
sub apply_pulse
{
	my ($self, $sigma, $spin, $hamil, $w1, $phase, $t) = @_;
	my ($res, $spinn, $sp, $sp_sys);
	my ($matx, $maty, $mat, $ham, $hamilt);
	my @ap_spin;

	$spinn = scalar(@{$self->{spin_array}});
	local $_ = $spin; @ap_spin = split;

	if (scalar(@ap_spin) > $spinn)
	{
		croak "The number of spins is inconsistent with the spin system [apply_pulse]";
	}

	$ham = undef;
	foreach $sp (@ap_spin)
	{
		$matx = "";
		$maty = "";

		foreach $sp_sys (@{$self->{spin_array}})
		{
			if ($sp eq $sp_sys)
			{
				$matx = "$matx" . "x ";
				$maty = "$maty" . "y ";
			}
			else
			{
				$matx = "$matx" . "e ";
				$maty = "$maty" . "e ";
			}
		}

		unless ($matx =~ /x/)
		{
			croak "The spin $sp does not exist in the spin system [apply_pulse]";
		}

		chop($matx); chop($maty);
		$mat = $self->make_density_matrix($w1*cos($phase), $matx, $w1*sin($phase), $maty);

		if (defined $ham)
		{
			$ham = $ham->add($mat);
		}
		else
		{
			$ham = $mat;
		}
	}

	if ($hamil == 0)
	{
		$hamilt = $ham->copy_matrix;
	}
	else
	{
		$hamilt = $ham->add($hamil); 
	}
	$res = $self->evolve_density_matrix($sigma, $hamilt, $t);
	return $res;
}

#U $res = $spin_sys->apply_pulse2($sigma, "I", $hamil, 2*$PI*3000, $PI/2, 10.0e-6);
# apply pulses to two spins at the same time
# "I I I" is not allowd. It must be "H S K", for example.
sub apply_pulse2
{
	my ($self, $sigma, $spin, $hamil, $w1, $phase, $t) = @_;
	my ($res, $spinn, $sp, $sp_sys);
	my ($matx, $maty, $mat, $ham, $hamilt, $rp);
	my @ap_spin;

	$spinn = scalar(@{$self->{spin_array}});
	local $_ = $spin; @ap_spin = split;

	if (scalar(@ap_spin) > $spinn)
	{
		croak "The number of spins is inconsistent with the spin system [apply_pulse2]";
	}

	$ham = undef;
	foreach $sp (@ap_spin)
	{
		$matx = "";
		$maty = "";
		$rp = 0;

		foreach $sp_sys (@{$self->{spin_array}})
		{
			if (($sp eq $sp_sys) && ($rp == 0))
			{
				$matx = "$matx" . "x ";
				$maty = "$maty" . "y ";
				$rp++;
			}
			else
			{
				$matx = "$matx" . "e ";
				$maty = "$maty" . "e ";
			}
		}

		unless ($matx =~ /x/)
		{
			croak "The spin $sp does not exist in the spin system [apply_pulse]";
		}

		chop($matx); chop($maty);
		$mat = $self->make_density_matrix($w1*cos($phase), $matx, $w1*sin($phase), $maty);

		if (defined $ham)
		{
			$ham = $ham->add($mat);
		}
		else
		{
			$ham = $mat;
		}
	}

	if ($hamil == 0)
	{
		$hamilt = $ham->copy_matrix;
	}
	else
	{
		$hamilt = $ham->add($hamil); 
	}
	$res = $self->evolve_density_matrix($sigma, $hamilt, $t);
	return $res;
}

1;
