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

#U $net = TI_AI->new(1, 784, 50, 100, 10);
# $p = TI_AI->new(2, 3, 4);
# ${$p->{wm}}[0]->show_matrix("w0: ");
# print("w0[0, 0]imaginary: ${$p->{wm}}[0]->[0][0][1]\n");
# $p->{inp}->show_matrix("inp: ");
sub new
{
	my $class = shift;
	my @mats = @_;
	my $self = {};	# the reference to an anonymous hash
	my ($i, $j, $k, @zr, @zrr, $inp, $oup, $w, $b, $wm, $bm);

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

		if ($i == 1)
		{
			$inp = $w;
			# $inp->show_matrix("inp: ");
		}
		else
		{
			push(@{$wm}, $w);
			# $w->show_matrix("w$j: ");

			@zrr = ();
			for ($k = 0; $k < $mats[0]; $k++)
			{
				@zr = (0) x $mats[$i];
				push(@zrr, [@zr]);
			}
			$b = TI_matrix->new(@zrr);
			push(@{$bm}, $b);
			# $b->show_matrix("b$j: ");
		}
	}

	$oup = TI_matrix->new(@zrr);
	# $oup->show_matrix("oup: ");

	$self->{wm} = $wm;
	$self->{bm} = $bm;
	$self->{inp} = $inp;
	$self->{oup} = $oup;

	bless($self, $class);
	return $self;
}

#U $y = $net->forward_dot;
# (a1, a2, a3) = (x1, x2).({w11, w12, w13}, {w21, w22, w23}) + (b1, b2, b3)
# not (a1, a2, a3) = (1, x1, x2).({b1, b2, b3}, {w11, w12, w13}, {w21, w22, w23})
sub forward_dot
{
	my $p = shift;
	my $m;

	$m = TI_matrix->new([1.0, 0.5]);
	$p->{inp}->mcopy($m);

	$m = TI_matrix->new([0.1, 0.3, 0.5], [0.2, 0.4, 0.6]);
	${$p->{wm}}[0]->mcopy($m);

	$m = TI_matrix->new([0.1, 0.2, 0.3]);
	${$p->{bm}}[0]->mcopy($m);

	$m = TI_matrix->new([0.1, 0.4], [0.2, 0.5], [0.3, 0.6]);
	${$p->{wm}}[1]->mcopy($m);

	$m = TI_matrix->new([0.1, 0.2]);
	${$p->{bm}}[1]->mcopy($m);

	$m = TI_matrix->new([0.1, 0.3], [0.2, 0.4]);
	${$p->{wm}}[2]->mcopy($m);

	$m = TI_matrix->new([0.1, 0.2]);
	${$p->{bm}}[2]->mcopy($m);

	$m = $p->{inp}->multiply(${$p->{wm}}[0]);
	$m = $m->add(${$p->{bm}}[0]);
	$m = $m->sigmoid;

	$m = $m->multiply(${$p->{wm}}[1]);
	$m = $m->add(${$p->{bm}}[1]);
	$m = $m->sigmoid;

	$m = $m->multiply(${$p->{wm}}[2]);
	$m = $m->add(${$p->{bm}}[2]);

	return $m;
}

#U $w = $p->simple_conv_net;
# $p = TI_AI->new(1, 28*28, 30,,,);
# $w = $p->simple_conv_net;
sub simple_conv_net
{
	my $p = shift;
	my $m;

	# input: an image with a size of 28*28 = 784

	# convolution filter with a number of 30
	$m = TI_matrix->new([1.0, 0.5]);
	$p->{inp}->mcopy($m);

	return 1;	# ############################
}

1;
