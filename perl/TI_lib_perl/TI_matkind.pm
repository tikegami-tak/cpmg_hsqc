=head1 NAME

TI_matkind - Judge the kind of the matrix

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

use strict;
use Carp;

my $THRESH_Z = 1.0e-7;	# You can change this value.

#U $a->is_zero_matrix;
# if A is zero matrix, return 1. Otherwise, return 0.
sub is_zero_matrix
{
	my $self = shift;
	my ($row, $col);

	foreach $row (@{$self})
	{
		foreach $col (@{$row})
		{
			if ((abs($col->[0]) > $THRESH_Z) || (abs($col->[1]) > $THRESH_Z)) { return 0; }
		}
	}

	return 1;
}

#U $a->is_c_identity_matrix;
# if A is c*E matrix, return c. Otherwise, return 0.
sub is_c_identity_matrix
{
	my $self = shift;
	my ($row, $col, $nrow, $ncol);
	my ($re, $im, $abre, $abim, $tmp);
	my $retc = 0;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = $self->[$row][$col][0];
			$im = $self->[$row][$col][1];
			$abre = abs($re);
			$abim = abs($im);

			if ($abim > $THRESH_Z) { return 0; }	# Imaginary = 0

			if ($row == $col)
			{
				if ($row == 0) { $retc = $re; }
				else
				{
					$tmp = $re - $retc;
					$tmp = abs($tmp);
					if ($tmp > $THRESH_Z) { return 0; }
				}
			}
			else
			{
				if ($abre > $THRESH_Z) { return 0; }
			}
		}
	}

	return $retc;
}

#U $a->is_diagonal_matrix;
# if A is diagonal matrix, return 1. Otherwise, return 0.
sub is_diagonal_matrix
{
	my $self = shift;
	my ($row, $col, $nrow, $ncol);
	my ($re, $im, $abre, $abim);

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		for ($col = 0; $col < $ncol; $col++)
		{
			if ($row != $col)
			{
				$re = $self->[$row][$col][0];
				$im = $self->[$row][$col][1];
				$abre = abs($re);
				$abim = abs($im);

				if (($abre > $THRESH_Z) || ($abim > $THRESH_Z)) { return 0; }
			}
		}
	}

	return 1;
}

#U $a->is_hermitian_matrix;
# if A is hermitian, return 1. Otherwise, return 0.
sub is_hermitian_matrix
{
	my $self = shift;
	my ($herm, $dif);

	$herm = $self->hermitian_conjugate;
	$dif = $herm->subtract($self);
	if ($dif->is_zero_matrix) { return 1; }
	else { return 0; }
}

1;
