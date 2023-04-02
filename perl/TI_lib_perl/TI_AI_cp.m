=head1 NAME

TI_AI - AI main function

=head1 AUTHOR

Takahisa Ikegami <ikegamit@yokohama-cu.ac.jp>

=cut

package TI_AI;

use TI_matrix;
# use strict;
use Carp;
use Math::Complex;

#U $net = TI_AI->new(784, 50, 100, 10);
# $p = TI_AI->new(2, 3, 4);
# ${$p->{wm}}[0]->show_matrix("w0: ");
# print("w0[0, 0]real: ${$p->{wm}}[0]->[0][0][0]\n");
# $p->{inp}->show_matrix("inp: ");
sub new
{
	my $class = shift;
	my @mats = @_;
	my $self = {};	# the reference to an anonymous hash
	my ($i, $j, $k, @zr, @zrr, $inp, $w, $b, $wm, $bm);

	@zr = (0) x $mats[0];	# one dimension
	$inp = TI_matrix->new([@zr]);
	# $inp->show_matrix("inp: ");

	$wm = [];	# the reference to an anonymous array
	$bm = [];
	for ($i = 1; $i < scalar(@mats); $i++)
	{
		@zrr = ();
		$j = $i - 1;
		for ($k = 0; $k < $mats[$j]; $k++)
		{
			@zr = (0) x $mats[$i];
			push(@zrr, [@zr]);
		}
		$w = TI_matrix->new(@zrr);
		# $w->show_matrix("w$i: ");
		push(@{$wm}, $w);

		@zr = (0) x $mats[$i];
		$b = TI_matrix->new([@zr]);
		# $b->show_matrix("b$i: ");
		push(@{$bm}, $b);
	}

	$j = $i - 1;
	@zr = (0) x $mats[$j];
	$oup = TI_matrix->new([@zr]);
	# $oup->show_matrix("oup: ");

	$self->{wm} = $wm;
	$self->{bm} = $bm;
	$self->{inp} = $inp;
	$self->{oup} = $oup;

	bless($self, $class);
	return $self;
}

#U $y = $net->forward_dot($x, $w, $b);
# (a1, a2, a3) = (x1, x2).({w11, w12, w13}, {w21, w22, w23}) + (b1, b2, b3)
# not (a1, a2, a3) = (1, x1, x2).({b1, b2, b3}, {w11, w12, w13}, {w21, w22, w23})
sub forward_dot
{
	my $self = shift;
	my $x = shift;
	my $w = shift;
	my $b = shift;
	my ($c, $d);

$p = TI_AI->new(2, 3, 4);
$w0 = $p->{inp}->multiply(${$p->{wm}}[0]);
$w1 = $w0->add(${$p->{bm}}[0]);
$w1->show_matrix("w: ");

	return ($d);
}

1;
