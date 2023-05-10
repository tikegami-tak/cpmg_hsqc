#!/usr/bin/perl

# Process all the SER files at once.
# Takahisa Ikegami on 040631.

for ($i = 1; $i <= 18; $i++)
{
	$inpfile = "ser_" . "$i";
	$oupfile = "Cpmg_" . "$i" . ".ft";

	system("convft.com ../$inpfile $oupfile");
}
