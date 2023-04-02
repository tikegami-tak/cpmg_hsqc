#!/usr/bin/perl

print "Input search word:";
$srch = <STDIN>;
chomp $srch;

opendir(DIR, ".");
@file = readdir(DIR);
close DIR;

foreach $ef (@file)
{
	open(FD, "<$ef");
	next if -d FD;
	next if -B FD;
	$num = 0;
	while (<FD>)
	{
		chomp;
		$line = $_;
		if (/$srch/)
		{
			@word = split;
			if ($word[0] ne '//')
			{
				$num++;
			}
		}
	}
	if ($num > 0) { print "$ef\t\t$num times\n"; }
	close FD;
}
