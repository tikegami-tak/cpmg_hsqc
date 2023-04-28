#!/usr/bin/perl

# Separate the ser file into n parts.
# Takahisa Ikegami on 040615.

# --------------------------------------------------

$num_dat = 18;

$inpfile = "ser";

$eachbyte = 2048;               # how many TD number?

# --------------------------------------------------

open(INP, "< $inpfile") || die "Can not open the input file for read.\n";

for ($i = 1; $i <= $num_dat; $i++)
{
	my $fh = "OUTP${i}";
	open(${$fh}, "> ser_${i}")      || die "Can not open the output file $i for write.\n";
}

$cnt = 1;
$cntbyte = 0;
while (1)
{
        $rebyte = sysread(INP, $iidata, 4);
        if (!defined $rebyte)
        {
                die " Can not read the input data ## error.\n";
        }
        elsif ($rebyte == 0)
        {
                last;
        }
        elsif ($rebyte != 4)
        {
                die " Can not read 4 bytes ## error.\n";
        }

        $sem = $cnt % $num_dat;
	if ($sem == 0) { $sem = $num_dat; }

	for ($i = 1; $i <= $num_dat; $i++)
	{
        	if ($sem == $i)
		{
			my $fh = "OUTP${i}";
			$wrbyte = syswrite(${$fh}, $iidata, 4);
        		if (!defined $wrbyte)
        		{
                		die " Can not write the data ## error.\n";
        		}
		}
	}

        $cntbyte++;
        if ($cntbyte == $eachbyte)
        {
                $cnt++;
                $cntbyte = 0;
        }
}

close INP;
for ($i = 1; $i <= $num_dat; $i++)
{
	my $fh = "OUTP${i}";
	close(${$fh});
}
