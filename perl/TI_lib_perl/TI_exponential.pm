=head1 NAME

TI_exponential - Exponential matrix calculations

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

use strict;
use Carp;

my $THRESH_Z = 1.0e-7;	# You can change this value.
my $DEBUG = 0;			# debug mode

#U $b = $a->exponential_operator;
# B = exp(A)
sub exponential_operator
{
	my $self = shift;
	my ($nrow, $ncol, $cons, $retc);
	my ($row, $col, $re, $im);
	my ($emat, $res, $tmp, $tmp2, $aa);
	my ($rtmp, $ctmp);

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if ($nrow != $ncol)
	{
		croak "The row and column sizes must be the same [exponential_operator]";
	}

	$res = $self->copy_matrix;				# A
	$retc = $res->is_diagonal_matrix;		# A = diagonal ?

	if ($retc)
	{
		for ($row = 0; $row < $nrow; $row++)
		{
			for ($col = 0; $col < $ncol; $col++)
			{
				if ($row == $col)
				{
					$re = $self->[$row][$col][0];
					$im = $self->[$row][$col][1];
					$res->[$row][$col][0] = exp($re) * cos($im);
					$res->[$row][$col][1] = exp($re) * sin($im);
				}
				else
				{
					$res->[$row][$col][0] = 0.0;
					$res->[$row][$col][1] = 0.0;
				}
			}
		}
		if ($DEBUG) { print STDERR "EXPonential -- diagonal --\n"; }
	}
	else		# A is not a diagonal matrix.
	{
		$aa = $res->multiply($self);			# AA
		$retc = $aa->is_c_identity_matrix;		# AA = cE ?

		if ($retc < -"$THRESH_Z")					# AA = cE && c < 0
		{
			$emat = $self->make_identity_matrix;	# E
			# E cos(sqrt[-c]) + A sin(sqrt[-c]) / sqrt[-c]
			$rtmp = sqrt(-"$retc");
			$ctmp = cos($rtmp);
			$cons = sin($rtmp) / $rtmp;
			$tmp = $self->scalar_multiply($cons);
			$tmp2 = $emat->scalar_multiply($ctmp);
			$res = $tmp2->add($tmp);
			if ($DEBUG) { print STDERR "EXPonential -- Taylor cos-sin --\n"; }
		}
		else	# Taylor series
		{
			# $res = $self->exp_taylor_expansion;
			$res = $self->exp_double_angle;
			# $res->show_matrix("Taylor exp ...", 10);
		}
	}

	return $res;
}

#U $b = $a->exp_taylor_expansion;
# B = exp(A) by Taylor series
sub exp_taylor_expansion
{
	my $self = shift;
	my ($nrow, $ncol, $i, $cons, $acum);
	my ($res, $aaa, $tmp);
	my $order = 500;			# You can change the order.

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if ($nrow != $ncol)
	{
		croak "The row and column sizes must be the same [exp_taylor_expansion]";
	}

	$res = $self->make_identity_matrix;		# E
	$acum = 1;

	foreach $i (1..$order)
	{
		$acum *= $i;
		$cons = 1.0 / $acum;
		if ($i == 1) { $aaa = $self->copy_matrix; }
		else { $aaa = $aaa->multiply($self); }
		$tmp = $aaa->scalar_multiply($cons);
		$res = $res->add($tmp);
		# $res->show_matrix(" RES is $i\n");

		if ($DEBUG)
		{
			printf STDERR "  RES order is %d, norm is %f\n", $i, $tmp->infinity_norm;
		}

		if ($tmp->is_zero_matrix)
		{
			if ($DEBUG) { print STDERR "EXPonential -- Taylor expansion order $i --\n"; }
			last;
		}
		if ($i == $order)
		{
			print STDERR "EXPonential -- Taylor expansion order $i last --\n";
			croak "The exponential did not converge [exp_taylor_expansion]";
		}
	}

	return $res;
}

#U $b = $a->exp_double_angle;
# B = exp(A) by repeated application of the double angle formulae
sub exp_double_angle
{
	my $self = shift;
	my ($nrow, $ncol, $norm, $i, $k, $dang, $cons);
	my ($res, $aaa, $tmp, $amat, $cosm, $sinm, $emat, $acum);
	my $order = 500;			# You can change the order.

	$nrow = scalar(@{$self});
	$ncol = scalar(@{$self->[0]});

	if ($nrow != $ncol)
	{
		croak "The row and column sizes must be the same [exp_double_angle]";
	}

	$amat = $self->scalar_multiply("0 -1");	# exp(A)=cos(-iA)+i*sin(-iA)
	$norm = $amat->infinity_norm;

	$k = 0;	$dang = 1;
	while ($norm > $dang)
	{
		$dang *= 2;
		$k++;
	}
	if ($k == 0)
	{
		$res = $self->exp_taylor_expansion;
		return $res;
		if ($DEBUG) { print STDERR "|A| ~= 1\n"; }
	}

	# Infinity norm of (-iA) is about 2^k
	if ($DEBUG) { print STDERR "|A| ~= 2^$k\n"; }

	$cons = 1.0 / $dang;
	$amat = $amat->scalar_multiply($cons);
	$cosm = $amat->make_identity_matrix;	# E
	$emat = $cosm->copy_matrix;
	$acum = 1;

	foreach $i (1..$order)
	{
		$acum *= $i;
		$cons = 1.0 / $acum;
		if ($i == 1) { $aaa = $amat->copy_matrix; }
		else { $aaa = $aaa->multiply($amat); }

		$tmp = $aaa->scalar_multiply($cons);

		if ($i % 2)	# for sin
		{
			if ($i == 1) { $sinm = $tmp->copy_matrix; }
			else
			{
				if (($i % 4) == 1) { $sinm = $sinm->add($tmp); }
				else { $sinm = $sinm->subtract($tmp); }
			}
		}
		else	# for cos
		{
			if ($i % 4) { $cosm = $cosm->subtract($tmp); }
			else { $cosm = $cosm->add($tmp); }
		}

		# $res->show_matrix(" RES is $i\n");

		if ($tmp->is_zero_matrix)
		{
			if ($DEBUG) { print STDERR "EXPonential -- Taylor double angle order $i --\n"; }
			last;
		}
		if ($i == $order)
		{
			print STDERR "EXPonential -- Taylor double angle order $i last --\n";
			croak "The exponential did not converge [exp_double_angle]";
		}
	}

	foreach $i (1..$k)
	{
		$tmp = $sinm->multiply($cosm);
		$sinm = $tmp->scalar_multiply(2);	# sin(2A)=2sin(A)cos(A)
		$tmp = $cosm->multiply($cosm);
		$tmp = $tmp->scalar_multiply(2);
		$cosm = $tmp->subtract($emat);		# cos(2A)=2cos^2(A)-E
	}

	$tmp = $sinm->scalar_multiply("0 1");	# exp(A)=cos(-iA)+i*sin(-iA)
	$res = $cosm->add($tmp);

	return $res;
}

1;
