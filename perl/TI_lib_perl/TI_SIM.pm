=head1 NAME

TI_SIM - Ti_simulator main function

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

package TI_SIM;

use TI_spinsys;
use TI_matrix;
use TI_pulse;
use TI_vectmodel;
use TI_relax;
use strict;
use Carp;
use Math::Complex;

#U $spin_sys = TI_SIM->new("I S T");
sub new
{
	my $class = shift;
	my $spin = shift;
	my $self = {};
	my ($spin_num, @spin_array);

	local $_ = $spin; @spin_array = split;
	$self->{spin_array} = [@spin_array];		# Spin system
	$spin_num = scalar(@spin_array);

	if ($spin_num < 1)
	{
		croak "The number of nuclei must be more than 1 [new]";
	}

	bless($self, $class);
	return $self;
}

1;
