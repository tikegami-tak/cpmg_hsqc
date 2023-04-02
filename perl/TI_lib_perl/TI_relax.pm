=head1 NAME

TI_relax - Calculate relaxation matrix

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

use strict;
use Carp;

#U $spin->cal_j_relax($sigma, $obs, 'dipolar_hetero', 'I S', 'dip', 'csa', 'I', 'cs');
sub cal_j_relax
{
	my $self = shift;
	my $sigma = shift;
	my $obs = shift;
	my @word = @_;
	my @rel_hamil = ();
	my @rel_cons = ();
	my @rel_omega = ();
	my @rel_quanta = ();
	my @rel_lab = ();
	my @j_spect = ();
	my @j_omega = ();
	my @j_cons = ();
	my @j_lab = ();
	my ($lab, $h, $i, $j, $n, $k, $m, $re, @ap_spin, @tmp_ap_spin, $swt);
	my ($reu, $imu, $re_norm, $im_norm, $dcom, $rcons, $rlab);

	$i = 0;

	if (scalar(@word) % 3)
	{
		croak "The number of the arguments is wrong. [cal_j_relax]";
	}

	for ($h = 0; $h < scalar(@word); $h+=3)
	{
		if ($word[$h] eq 'dipol')
		{
			@ap_spin = ();
			($swt, $word[$h+1]) = $self->arrage_spin_lab($word[$h+1]);
			local $_ = $word[$h+2]; @tmp_ap_spin = split;
			if ($swt)
			{
				$ap_spin[0] = $tmp_ap_spin[1];
				$ap_spin[1] = $tmp_ap_spin[0];
			}
			else
			{
				$ap_spin[0] = $tmp_ap_spin[0];
				$ap_spin[1] = $tmp_ap_spin[1];
			}
	
			$lab = $self->expand_dens_mat_label($word[$h+1], 'z z');
			$rel_hamil[$i] = $self->make_density_matrix(2/sqrt(6), $lab);
			if ($ap_spin[0] ge $ap_spin[1])
			{
				$rel_cons[$i] = "D[$ap_spin[0],$ap_spin[1]]";
			}
			else
			{
				$rel_cons[$i] = "D[$ap_spin[1],$ap_spin[0]]";
			}
			$rel_lab[$i] = 'd';
			$rel_omega[$i] = '0';
			$rel_quanta[$i] = 0;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'p n');
			$rel_hamil[$i] = $self->make_density_matrix(-1/2/sqrt(6), $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			if ($ap_spin[0] eq $ap_spin[1])
			{
				$rel_omega[$i] = '0';
			}
			elsif ($ap_spin[0] gt $ap_spin[1])
			{
				$rel_omega[$i] = "$ap_spin[0]" . '-' . "$ap_spin[1]";
			}
			else
			{
				$rel_omega[$i] = "$ap_spin[1]" . '-' . "$ap_spin[0]";
			}
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = 0;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'n p');
			$rel_hamil[$i] = $self->make_density_matrix(-1/2/sqrt(6), $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = $rel_omega[$i-1];
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = 0;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'z p');
			$rel_hamil[$i] = $self->make_density_matrix(-1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = "$ap_spin[1]";
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = 1;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'z n');
			$rel_hamil[$i] = $self->make_density_matrix(1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = "$ap_spin[1]";
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = -1;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'p z');
			$rel_hamil[$i] = $self->make_density_matrix(-1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = "$ap_spin[0]";
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = 1;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'n z');
			$rel_hamil[$i] = $self->make_density_matrix(1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = "$ap_spin[0]";
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = -1;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'p p');
			$rel_hamil[$i] = $self->make_density_matrix(1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_lab[$i] = 'd';
			if ($ap_spin[0] eq $ap_spin[1])
			{
				$rel_omega[$i] = '2' . "$ap_spin[0]";
			}
			elsif ($ap_spin[0] gt $ap_spin[1])
			{
				$rel_omega[$i] = "$ap_spin[0]" . '+' . "$ap_spin[1]";
			}
			else
			{
				$rel_omega[$i] = "$ap_spin[1]" . '+' . "$ap_spin[0]";
			}
			$rel_quanta[$i] = 2;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'n n');
			$rel_hamil[$i] = $self->make_density_matrix(1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_omega[$i] = $rel_omega[$i-1];
			$rel_lab[$i] = 'd';
			$rel_quanta[$i] = -2;
			$i++;
		}
		elsif ($word[$h] eq 'csa')
		{
			local $_ = $word[$h+2]; @ap_spin = split;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'z');
			$rel_hamil[$i] = $self->make_density_matrix(2/sqrt(6), $lab);
			$rel_cons[$i] = "C[$ap_spin[0]]";
			$rel_lab[$i] = 'c';
			$rel_omega[$i] = '0';
			$rel_quanta[$i] = 0;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'p');
			$rel_hamil[$i] = $self->make_density_matrix(-1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_lab[$i] = 'c';
			$rel_omega[$i] = "$ap_spin[0]";
			$rel_quanta[$i] = 1;
			$i++;

			$lab = $self->expand_dens_mat_label($word[$h+1], 'n');
			$rel_hamil[$i] = $self->make_density_matrix(1/2, $lab);
			$rel_cons[$i] = $rel_cons[$i-1];
			$rel_lab[$i] = 'c';
			$rel_omega[$i] = "$ap_spin[0]";
			$rel_quanta[$i] = -1;
			$i++;
		}
		elsif ($word[$h] eq 'quad')
		{
			print "\n Q relaxation\n";
		}
		else
		{
			croak "The relax hamiltonian is not yet implemented. $word[$h] [make_relax_hamil]";
		}
	}

	$n = $i;
	$k = 0;

	for ($i = 0; $i < $n; $i++)
	{
		for ($j = 0; $j < $n; $j++)
		{
			if (($rel_omega[$i] eq $rel_omega[$j]) && ($rel_quanta[$i] == -$rel_quanta[$j]))
			{
				# <obs | [el_Apq, [Apq, sigma]] / <obs|obs>
				$dcom = $sigma->double_commutator($rel_hamil[$i], $rel_hamil[$j]);
				($reu, $imu, $re_norm, $im_norm) = $self->vector_projection($obs, $dcom);
				$re = $re_norm * 0.5;
				if ($rel_quanta[$i]%2) { $re = -$re; }
				if ($rel_cons[$i] gt $rel_cons[$j]) { $rcons = "$rel_cons[$i]" . " $rel_cons[$j]"; }
				else { $rcons = "$rel_cons[$j]" . " $rel_cons[$i]"; }
				if ($rel_lab[$i] gt $rel_lab[$j]) { $rlab = "$rel_lab[$i]" . "$rel_lab[$j]"; }
				else { $rlab = "$rel_lab[$j]" . "$rel_lab[$i]"; }

				if ($re)
				{
					$m = 0;
					for ($m = 0; $m < $k; $m++)
					{
						if (($j_omega[$m] eq $rel_omega[$i]) && ($j_cons[$m] eq $rcons))
						{
							$j_spect[$m] += $re;
							last;
						}
					}
					if ($m == $k)
					{
						$j_spect[$m] = $re;
						$j_omega[$m] = $rel_omega[$i];
						$j_cons[$m] = $rcons;
						$j_lab[$m] = $rlab;
						$k++;
					}
				}
			}
		}
	}

	print "\n";
	print 'd = u/4pi h/2pi ri rs / r^3';
	print "\n";
	print 'c = (Sh - Sv) ri B0 / root3';
	print "\n\n";

	for ($m = 0; $m < $k; $m++)
	{

		if ($j_lab[$m] eq 'dd')
		{
			printf "  J ( %s ) ... %s /8 * %.3f\n", $j_omega[$m], $j_cons[$m], $j_spect[$m]*48;
		}
		elsif ($j_lab[$m] eq 'cc')
		{
			printf "  J ( %s ) ... %s /6 * %.3f\n", $j_omega[$m], $j_cons[$m], $j_spect[$m]*12;
		}
		elsif (($j_lab[$m] eq 'cd') || ($j_lab[$m] eq 'dc'))
		{
			printf "  J ( %s ) ... %s /(4*root3) * %.3f\n", $j_omega[$m], $j_cons[$m], $j_spect[$m]*24;
		}
		else
		{
			croak "The constant is wrong. [make_relax_hamil]";
		}
	}
}

#U $spin_sys->arrage_spin_lab('N H');
# if the spin_sys is "I H N S", then "H N" will be returned.
sub arrage_spin_lab
{
	my ($self, $spin) = @_;
	my (@ap_spin, $sp, $i, $j, $n0, $n1);

	local $_ = $spin; @ap_spin = split;
	if (scalar(@ap_spin) != 2)
	{
		croak "The number of spins must be 2. [arrage_spin_lab]";
	}

	$n0 = undef; $n1 = undef;

	for ($j = 0; $j < 2; $j++)
	{
		$sp = $ap_spin[$j];

		for ($i = 0; $i < scalar(@{$self->{spin_array}}); $i++)
		{
#			if ($sp eq ${$self->{spin_array}}[$i])
			if ($sp eq $self->{spin_array}->[$i])
			{
				if ($j == 0)	{ $n0 = $i; last; }
				else 			{ $n1 = $i; last; }
			}
		}
		if ($i == scalar(@{$self->{spin_array}}))
		{
			croak "Can not find the spin label $sp. [arrage_spin_lab]";
		}
	}

	if ($n0 < $n1)
	{
		return (0, "$ap_spin[0] $ap_spin[1]");
	}
	elsif ($n1 < $n0)
	{
		return (1, "$ap_spin[1] $ap_spin[0]");
	}
	else
	{
		croak "The two spin labels are the same [arrage_spin_lab]";
	}
}

#U $spin_sys->expand_dens_mat_label("H N", "z x");
# if the spin_sys is "I H N S", then "e z x e" will be returned.
sub expand_dens_mat_label
{
	my $self = shift;
	my ($spin, $pd) = @_;
	my (@ap_spin, @ap_pd);
	my ($spinn, $mat, $sp_sys, $i);

	local $_ = $spin; @ap_spin = split;
	$spinn = scalar(@ap_spin);

	local $_ = $pd; @ap_pd = split;

	if ($spinn != scalar(@ap_pd))
	{
		croak "The numbers of spins in the argument are not the same [expand_dens_mat_label]";
	}
	if ($spinn > scalar(@{$self->{spin_array}}))
	{
		croak "The number of spins is inconsistent with the spin system [expand_dens_mat_label]";
	}

	$i = 0;
	$mat = "";

	foreach $sp_sys (@{$self->{spin_array}})
	{
		if ($i < $spinn)
		{
			if ($ap_spin[$i] eq $sp_sys)
			{
				$mat = "$mat" . "$ap_pd[$i] ";
				$i++;
				next;
			}
		}
		$mat = "$mat" . "e ";
	}
	chop($mat);

	return $mat;
}

1;
