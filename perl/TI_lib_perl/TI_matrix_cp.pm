=head1 NAME

TI_matrix - Matrix calculations

=head1 AUTHOR

Takahisa Ikegami <ikegamit@yokohama-cu.ac.jp>

=cut

package TI_matrix;

use TI_exponential;
use TI_matkind;
use strict;
use Carp;
use Math::Complex;

my $THRESH_Z = 1.0e-7;	# Can be changed.

#U $a = new TI_matrix ([1, 0, 2-4*i, -$PI], [rand, 3+4, 4*2-5, 0]);
sub new
{
	my $class = shift;
	my $self;
	my $ncol = scalar(@{$_[0]});
	my ($re, $im, $col, $row, $line, $val);

	$self = [[[]]];
	$row = 0;
	foreach $line (@_)		# $line is the reference to each row.
	{
		$self->[$row] = [[]];
		if (scalar(@{$line}) != $ncol)
		{
			croak "The column sizes are not the same [new]";
		}
		$col = 0;
		foreach $val (@{$line})
		{
			$self->[$row][$col] = [];

			if ($val =~ /i/)
			{
				$self->[$row][$col][0] = ($val->Re);
				$self->[$row][$col][1] = ($val->Im);
			}
			else
			{
				local $_ = $val;	($re, $im) = split;
				$im = 0 unless defined $im;
				$self->[$row][$col][0] = $re;	# Real
				$self->[$row][$col][1] = $im;	# Imaginary
			}

			$col++;
		}
		$row++;
	}

	bless $self, $class;
	return $self;
}

#U $c = $a->add($b);
# C = A + B
sub add
{
	my $self = shift;
	my $other = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if (($nrow != scalar(@{$other})) || ($ncol != scalar(@{$other->[0]})))
	{
		croak "The matrix sizes are inconsistent [add]";
	}

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = $self->[$row][$col][0] + $other->[$row][$col][0];
			$im = $self->[$row][$col][1] + $other->[$row][$col][1];
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $c = $a->subtract($b);
# C = A - B
sub subtract
{
	my $self = shift;
	my $other = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if (($nrow != scalar(@{$other})) || ($ncol != scalar(@{$other->[0]})))
	{
		croak "The matrix sizes are inconsistent [subtract]";
	}

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = $self->[$row][$col][0] - $other->[$row][$col][0];
			$im = $self->[$row][$col][1] - $other->[$row][$col][1];
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $c = $a->multiply($b);
# C = A * B
sub multiply
{
	my $self = shift;
	my $other = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im, $i);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$other->[0]});

	if (scalar(@{$self->[0]}) != scalar(@{$other}))
	{
		croak "The matrix sizes are inconsistent [multiply]";
	}

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = 0; $im = 0;
			for ($i = 0; $i < scalar(@{$other}); $i++)
			{
				$re += ($self->[$row][$i][0] * $other->[$i][$col][0]
						- $self->[$row][$i][1] * $other->[$i][$col][1]);

				$im += ($self->[$row][$i][0] * $other->[$i][$col][1]
						+ $self->[$row][$i][1] * $other->[$i][$col][0]);
			}
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $b = $a->scalar_multiply(3*(1+3*i));
# B = (3+9I) * A
sub scalar_multiply
{
	my $self = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my ($scalar_re, $scalar_im);
	my $resrow;

	if ($_[0] =~ /i/)
	{
		$scalar_re = ($_[0]->Re);
		$scalar_im = ($_[0]->Im);
	}
	else
	{
		local $_ = $_[0];	($scalar_re, $scalar_im) = split;
		$scalar_im = 0 unless defined $scalar_im;
	}

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = ($self->[$row][$col][0]) * $scalar_re
				- ($self->[$row][$col][1]) * $scalar_im;
			$im = ($self->[$row][$col][0]) * $scalar_im
				+ ($self->[$row][$col][1]) * $scalar_re;

			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $b = $a->hermitian_conjugate;
# B = hermitian conjugate of A
sub hermitian_conjugate
{
	my $self = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $ncol; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $nrow; $col++)
		{
			$re = $self->[$col][$row][0];
			$im = -"$self->[$col][$row][1]";
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $c = $a->complex_conjugate;
# C = complex conjugate of A
sub complex_conjugate
{
	my $self = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = $self->[$row][$col][0];
			$im = -"$self->[$row][$col][1]";
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $c = $a->transpose;
# C = transpose of A
sub transpose
{
	my $self = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $ncol; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $nrow; $col++)
		{
			$re = $self->[$col][$row][0];
			$im = $self->[$col][$row][1];
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U ($re, $im) = $a->trace;
# trace of A
sub trace
{
	my $self = shift;
	my ($nrow, $ncol, $i, $res_re, $res_im);

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if ($ncol != $nrow)
	{
		croak "The column and row sizes are not the same [trace]";
	}

	$res_re = 0; $res_im = 0;
	for ($i = 0; $i < $nrow; $i++)
	{
		$res_re += $self->[$i][$i][0];
		$res_im += $self->[$i][$i][1];
	}

	return ($res_re, $res_im);
}

#U $a->show_matrix("Any comment "); %.2f output.
#U $a->show_matrix("Any comment ", 4); %.4f output.
sub show_matrix
{
	my $self = shift;
	my ($row, $col);
	my ($nraw, $ncol);
	my ($mes, $od, $sign);

	$od = 2;
	if (@_ >= 2)
	{
		($mes, $od) = @_;
		print "$mes";
	}
	elsif (@_ == 1)
	{
		($mes) = @_;
		print "$mes";
	}

	$nraw = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	print " $nraw raw X $ncol column matrix\n\n";

	foreach $row (@{$self})
	{
		foreach $col (@{$row})
		{
			if (($col->[1]) < 0) { $sign = ' '; }
			else { $sign = ' +'; }
			printf("%.${od}f%s%.${od}f I	", $col->[0], $sign, $col->[1]);
		}
		print "\n";
	}
	print "\n";
}

#U $c = $a->direct_product($b);
# C = A (X) B
sub direct_product
{
	my $mata = shift;
	my $matb = shift;
	my @zres = ();
	my ($nrowa, $ncola, $nrowb, $ncolb);
	my ($row, $col);
	my ($re, $im, $i, $j, $k, $m);
	my ($resrow, $newres);

	$nrowa = scalar(@{$mata});
	$ncola = scalar(@{$mata->[0]});
	$nrowb = scalar(@{$matb});
	$ncolb = scalar(@{$matb->[0]});

	for ($row = 0; $row < ($nrowa * $nrowb); $row++)
	{
		$resrow = [];
		for ($col = 0; $col < ($ncola * $ncolb); $col++)
		{
			push(@{$resrow}, "0 0");
		}
		push(@zres, $resrow);
	}
	$newres = new TI_matrix @zres;

	for ($i = 0; $i < $nrowa; $i++)
	{
		for ($j = 0; $j < $ncola; $j++)
		{
			for ($k = 0; $k < $nrowb; $k++)
			{
				for ($m = 0; $m < $ncolb; $m++)
				{
					$re = ($mata->[$i][$j][0] * $matb->[$k][$m][0]
						- $mata->[$i][$j][1] * $matb->[$k][$m][1]);

					$im = ($mata->[$i][$j][0] * $matb->[$k][$m][1]
						+ $mata->[$i][$j][1] * $matb->[$k][$m][0]);

					$col = $ncolb * $j + $m;
					$row = $nrowb * $i + $k;

					$newres->[$row][$col][0] = $re;
					$newres->[$row][$col][1] = $im;
				}
			}
		}
	}

	return $newres;
}

#U $d = $a->double_commutator($b, $c);
# D = [C, [B, A]]
sub double_commutator
{
	my $oth_a = shift;
	my $oth_b = shift;
	my $oth_c = shift;
	my ($former, $latter);

	$former = $oth_b->commutator($oth_a);
	$latter = $oth_c->commutator($former);
	return ($latter);
}

#U $c = $a->commutator($b);
# C = [A, B]
sub commutator
{
	my $self = shift;
	my $other = shift;
	my ($nrow, $former, $latter);

	$nrow = scalar(@{$self});

	if (   ($nrow != scalar(@{$self->[0]}))
		|| ($nrow != scalar(@{$other}))
		|| ($nrow != scalar(@{$other->[0]})))
	{
		croak "The row and column sizes must be the same [commutator]";
	}

	$former = $self->multiply($other);
	$latter = $other->multiply($self);
	return ($former->subtract($latter));
}

#U $super_a = $a->commutator_superoperator;
# A (X) Et - E (X) At
sub commutator_superoperator
{
	my $self = shift;
	my ($emat, $temat, $tself, $tmp1, $tmp2, $nrow);

	$nrow = scalar(@{$self});

	if ($nrow != scalar(@{$self->[0]}))
	{
		croak "The row and column sizes must be the same [commutator_superoperator]";
	}

	$emat = $self->make_identity_matrix;
	$temat = $emat->transpose;		# the same as $emat
	$tself = $self->transpose;
	$tmp1 = $self->direct_product($temat);
	$tmp2 = $emat->direct_product($tself);
	return ($tmp1->subtract($tmp2));
}

#U $b = $a->copy_matrix;
# B = A
sub copy_matrix
{
	my $mat = shift;
	my ($nrow, $ncol, $row, $col);
	my ($resrow, $comp);
	my @zres = ();

	$nrow = scalar(@{$mat});
	$ncol = scalar(@{$mat->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$comp = "$mat->[$row][$col][0] $mat->[$row][$col][1]";
			push(@{$resrow}, $comp);
		}
		push(@zres, $resrow);
	}

	return (new TI_matrix @zres);
}

#U $e = $a->make_identity_matrix;
# make identity matrix E with the same size of A
sub make_identity_matrix
{
	my $self = shift;
	my @res = ();
	my $resrow;
	my ($nrow, $row, $col);

	$nrow = scalar(@{$self});

	if ($nrow != scalar(@{$self->[0]}))
	{
		croak "The row and column sizes must be the same [make_identity_matrix]";
	}

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $nrow; $col++)
		{
			if ($row == $col)
			{
				push(@{$resrow}, 1);
			}
			else
			{
				push(@{$resrow}, 0);
			}
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $b = $a->mat_to_column_vector;
# B is the column vector of A
sub mat_to_column_vector
{
	my $self = shift;
	my @res = ();
	my ($nrow, $ncol, $row, $col, $re, $im);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	for ($row = 0; $row < $nrow; $row++)
	{
		for ($col = 0; $col < $ncol; $col++)
		{
			$resrow = [];
			$re = $self->[$row][$col][0];
			$im = $self->[$row][$col][1];
			push(@{$resrow}, "$re $im");
			push(@res, $resrow);
		}
	}

	return (new TI_matrix @res);
}

#U $b = $a->column_vector_to_mat(3, 4);
# B is the (3, 4) matrix of a column vector A
sub column_vector_to_mat
{
	my $self = shift;
	my ($nrow, $ncol) = @_;
	my @res = ();
	my ($row, $col, $crow, $re, $im);
	my $resrow;

	if (scalar(@{$self->[0]}) != 1)
	{
		croak "The matrix must be a column vector [column_vector_to_mat]";
	}

	if (scalar(@{$self}) != ($nrow * $ncol))
	{
		croak "The column number times row number is inconsistent [column_vector_to_mat]";
	}

	for ($row = 0; $row < $nrow; $row++)
	{
		$resrow = [];
		for ($col = 0; $col < $ncol; $col++)
		{
			$crow = $row * $nrow + $col;
			$re = $self->[$crow][0][0];
			$im = $self->[$crow][0][1];
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}

	return (new TI_matrix @res);
}

#U $norm = $a->infinity_norm;
# Return infinity-norm of the matrix
sub infinity_norm
{
	my $self = shift;
	my ($row, $col, $nrow, $ncol, $re, $im);
	my ($sqn, $sum_sqn, $max_abs);

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	$max_abs = 0;
	for ($row = 0; $row < $nrow; $row++)
	{
		$sum_sqn = 0.0;
		for ($col = 0; $col < $ncol; $col++)
		{
			$re = $self->[$row][$col][0];
			$im = $self->[$row][$col][1];
			$sqn = sqrt($re * $re + $im * $im);
			$sum_sqn += $sqn;
		}

		if ($max_abs < $sum_sqn) { $max_abs = $sum_sqn; }
	}

	return $max_abs;
}

#U ($re, $im) = $a->show_matrix_ele(3, 4);
sub show_matrix_ele
{
	my $self = shift;
	my ($row, $col) = @_;
	my ($nraw, $ncol, $re, $im);

	$nraw = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if (($row >= $nraw) || ($col >= $ncol))
	{
		croak "The element is out of the matrix [show_matrix_ele]";
	}

	$re = $self->[$row][$col][0];
	$im = $self->[$row][$col][1];

	return ($re, $im);
}

#U $c = $a->convol($b, $step);
# C = A ** B(filter)
sub convol
{
	my $self = shift;
	my $filt = shift;
	my $s = shift;
	my @res = ();
	my ($prow, $pcol, $nrow, $ncol, $frow, $fcol, $orow, $ocol, $row, $col, $re, $im, $s, $r, $c, $gr, $gc, $rev, $imv);
	my $resrow;

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	$frow = scalar(@{$filt});
	$fcol = scalar(@{$filt->[0]});

	# if ((($nrow + 2*$p -$frow) % $s != 0) || (($ncol + 2*$p -$fcol) % $s != 0))
	# {
	#	croak "The matrix sizes are inconsistent [convol]";
	# }

	$prow = int($frow/2);
	$pcol = int($fcol/2);

	# If the filter size is odd, the size of input is the same as that of output.
	$orow = int(($nrow + 2*$prow -$frow)/$s +1.1);
	$ocol = int(($ncol + 2*$pcol -$fcol)/$s +1.1);

	for ($row = 0; $row < $orow; $row+=$s)
	{
		$resrow = [];
		for ($col = 0; $col < $ocol; $col+=$s)
		{
			$re = 0;
			$im = 0;
			for ($r = 0; $r < $frow; $r++)
			{
				$gr = $row - $prow + $r;
				for ($c = 0; $c < $fcol; $c++)
				{
					$gc = $col - $pcol + $c;

					if (($gr >= 0) && ($gr < $nrow) && ($gc >= 0) && ($gc < $ncol))
					{
						$rev = $self->[$gr][$gc][0];
						$imv = $self->[$gr][$gc][1];
					}
					else
					{
						$rev = 0;
						$imv = 0;
					}
					$re += ($rev * $filt->[$r][$c][0] - $imv * $filt->[$r][$c][1]);
					$im += ($rev * $filt->[$r][$c][1] + $imv * $filt->[$r][$c][0]);
				}
			}
			push(@{$resrow}, "$re $im");
		}
		push(@res, $resrow);
	}
	return (new TI_matrix @res);
}

1;
