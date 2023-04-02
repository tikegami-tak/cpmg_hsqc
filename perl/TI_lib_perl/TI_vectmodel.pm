=head1 NAME

TI_vectmodel - Draw vector model

=head1 AUTHOR

Takahisa Ikegami <ikegami@bs.aist-nara.ac.jp>

=cut

use strict;
use Carp;

#U ($theta, $phi) = $spin_sys->vector_model_xyz($sigma, "x e", "y e", "z e");
sub vector_model_xyz
{
	my ($self, $sigma, $x_dens, $y_dens, $z_dens) = @_;
	my ($obs, $dens, $xy, $theta, $phi);
	my ($reu, $imu, $re_norm, $im_norm);
	my @xyz = ();

	foreach $dens ($x_dens, $y_dens, $z_dens)
	{
		$obs = $self->make_density_matrix($dens);
		($reu, $imu, $re_norm, $im_norm) = $self->vector_projection($obs, $sigma);
		push(@xyz, $re_norm);	# I do not know where the imaginary part should go.
	}

	$xy = sqrt($xyz[0] * $xyz[0] + $xyz[1] * $xyz[1]);

	# angle between the x axis and the projection of the vector on the xy plane.
	$phi = atan2($xyz[1], $xyz[0]);
	# angle beteen the vector and the projection of the vector on the xy plane.
	$theta = atan2("$xyz[2]", "$xy");

	# If the double quotation marks are removed when $xy is euqal to 0,
	# atan2 does not seem to return the real answer.
	# This may be a bug in the perl mathematics library.

	return ($theta, $phi);
}

#U ($reu, $imu, $re_norm, $im_norm) = $spin_sys->vector_projection($obs, $sigma);
# <obs|sigma> / <obs|obs>
sub vector_projection
{
	my ($self, $obs, $sigma) = @_;
	my ($el_obs, $tmp);
	my ($reu, $imu, $red, $imd, $re_norm, $im_norm);

	$el_obs = $obs->hermitian_conjugate;
	$tmp = $el_obs->multiply($sigma);
	($reu, $imu) = $tmp->trace;
	$tmp = $el_obs->multiply($obs);
	($red, $imd) = $tmp->trace;

	# (reu + imu I) / (red + imd I)
	$re_norm = ($reu * $red + $imu * $imd) / ($red * $red + $imd * $imd);
	$im_norm = ($imu * $red - $reu * $imd) / ($red * $red + $imd * $imd);

	return ($reu, $imu, $re_norm, $im_norm);
}

1;
