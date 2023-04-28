#!/usr/bin/perl

#------------------------------------------------
$sta_res = 655-14;
$end_res = 699+100;
#------------------------------------------------

$line = 'MASMTGGQQMGRGSAWQVNTAYTAGQLVTYNGKTYKCLQPHTSLAGWEPSNVPALWQLQxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx';

@res_nam = split //, $line;

for ($i = $sta_res; $i <= $end_res; $i++)
{
	$j = $i - $sta_res;
	$amide[$j] = '?';
	$hamide[$j] = '?';
}

open(SHIFT, "< assign_trosy4.list") || die "Cannot open the Sparky assignment file.\n";

SLINE:while (<SHIFT>)
{
	$line = $_;
	next SLINE if (/^#/ || /^!/ || /^\n/);
	chomp;
	@word = split;

	# $word[8] =~ /^(\d+)(\D)-(\d+)(\D)$/;

	if ($word[0] =~ /Ne/)	# side-chain
	{
		$word[0] =~ /^(\D)(\d+)Ne-HNe$/;
		$res = $1;
		$j = $2 - $sta_res + 100;
	}
	else
	{
		$word[0] =~ /^(\D)(\d+)N-HN$/;
		$res = $1;
		$j = $2 - $sta_res;
	}
	$amide[$j] = $word[1];
	$hamide[$j] = $word[2];

	# if    ($word[2] =~ /^ALA/)	{ $res = 'A'; }
	# elsif ($word[2] =~ /^CYS/)	{ $res = 'C'; }
	# elsif ($word[2] =~ /^ASP/)	{ $res = 'D'; }
	# elsif ($word[2] =~ /^GLU/)	{ $res = 'E'; }
	# elsif ($word[2] =~ /^PHE/)	{ $res = 'F'; }
	# elsif ($word[2] =~ /^GLY/)	{ $res = 'G'; }
	# elsif ($word[2] =~ /^HIS/)	{ $res = 'H'; }
	# elsif ($word[2] =~ /^ILE/)	{ $res = 'I'; }
	# elsif ($word[2] =~ /^LYS/)	{ $res = 'K'; }
	# elsif ($word[2] =~ /^LEU/)	{ $res = 'L'; }
	# elsif ($word[2] =~ /^MET/)	{ $res = 'M'; }
	# elsif ($word[2] =~ /^ASN/)	{ $res = 'N'; }
	# elsif ($word[2] =~ /^PRO/)	{ $res = 'P'; }
	# elsif ($word[2] =~ /^GLN/)	{ $res = 'Q'; }
	# elsif ($word[2] =~ /^ARG/)	{ $res = 'R'; }
	# elsif ($word[2] =~ /^SER/)	{ $res = 'S'; }
	# elsif ($word[2] =~ /^THR/)	{ $res = 'T'; }
	# elsif ($word[2] =~ /^VAL/)	{ $res = 'V'; }
	# elsif ($word[2] =~ /^TRP/)	{ $res = 'W'; }
	# elsif ($word[2] =~ /^TYR/)	{ $res = 'Y'; }

	# if    ($word[0] =~ /^Ala/)	{ $res[$i] = 'A'; }
	# elsif ($word[0] =~ /^Cys/)	{ $res[$i] = 'C'; }
	# elsif ($word[0] =~ /^Asp/)	{ $res[$i] = 'D'; }
	# elsif ($word[0] =~ /^Glu/)	{ $res[$i] = 'E'; }
	# elsif ($word[0] =~ /^Phe/)	{ $res[$i] = 'F'; }
	# elsif ($word[0] =~ /^Gly/)	{ $res[$i] = 'G'; }
	# elsif ($word[0] =~ /^His/)	{ $res[$i] = 'H'; }
	# elsif ($word[0] =~ /^Ile/)	{ $res[$i] = 'I'; }
	# elsif ($word[0] =~ /^Lys/)	{ $res[$i] = 'K'; }
	# elsif ($word[0] =~ /^Leu/)	{ $res[$i] = 'L'; }
	# elsif ($word[0] =~ /^Met/)	{ $res[$i] = 'M'; }
	# elsif ($word[0] =~ /^Asn/)	{ $res[$i] = 'N'; }
	# elsif ($word[0] =~ /^Pro/)	{ $res[$i] = 'P'; }
	# elsif ($word[0] =~ /^Gln/)	{ $res[$i] = 'Q'; }
	# elsif ($word[0] =~ /^Arg/)	{ $res[$i] = 'R'; }
	# elsif ($word[0] =~ /^Ser/)	{ $res[$i] = 'S'; }
	# elsif ($word[0] =~ /^Thr/)	{ $res[$i] = 'T'; }
	# elsif ($word[0] =~ /^Val/)	{ $res[$i] = 'V'; }
	# elsif ($word[0] =~ /^Tyr/)	{ $res[$i] = 'Y'; }
	# elsif ($word[0] =~ /^Trps/)	{ $res[$i] = 'w'; }	# side-chain
	# elsif ($word[0] =~ /^Trp/)	{ $res[$i] = 'W'; }

	if    ($res =~ /^A/)	{ $res3nam[$j] = 'Ala'; }
	elsif ($res =~ /^C/)	{ $res3nam[$j] = 'Cys'; }
	elsif ($res =~ /^D/)	{ $res3nam[$j] = 'Asp'; }
	elsif ($res =~ /^E/)	{ $res3nam[$j] = 'Glu'; }
	elsif ($res =~ /^F/)	{ $res3nam[$j] = 'Phe'; }
	elsif ($res =~ /^G/)	{ $res3nam[$j] = 'Gly'; }
	elsif ($res =~ /^H/)	{ $res3nam[$j] = 'His'; }
	elsif ($res =~ /^I/)	{ $res3nam[$j] = 'Ile'; }
	elsif ($res =~ /^K/)	{ $res3nam[$j] = 'Lys'; }
	elsif ($res =~ /^L/)	{ $res3nam[$j] = 'Leu'; }
	elsif ($res =~ /^M/)	{ $res3nam[$j] = 'Met'; }
	elsif ($res =~ /^N/)	{ $res3nam[$j] = 'Asn'; }
	elsif ($res =~ /^P/)	{ $res3nam[$j] = 'Pro'; }
	elsif ($res =~ /^Q/)	{ $res3nam[$j] = 'Gln'; }
	elsif ($res =~ /^R/)	{ $res3nam[$j] = 'Arg'; }
	elsif ($res =~ /^S/)	{ $res3nam[$j] = 'Ser'; }
	elsif ($res =~ /^T/)	{ $res3nam[$j] = 'Thr'; }
	elsif ($res =~ /^V/)	{ $res3nam[$j] = 'Val'; }
	elsif ($res =~ /^W/)	{ $res3nam[$j] = 'Trp'; }
	elsif ($res =~ /^Y/)	{ $res3nam[$j] = 'Tyr'; }

	if ($res_nam[$j] ne $res)
	{
		$i = $j+$sta_res;
		print "Residue $i $res is different from $res_nam[$j].\n";
		$res3nam[$j] = 'Trp';
	}
}
close SHIFT;

open(OUP, "> ChiA1_shift.txt") || die "Cannot open the output file.\n";

for ($i = $sta_res; $i <= $end_res; $i++)
{
	$j = $i - $sta_res;

	if (($amide[$j] =~ /\w/) && ($hamide[$j] =~ /\w/))
	{
		printf OUP "%d   %s     %s	%.3f	%.3f\n",
			$i, $res3nam[$j], $res_nam[$j], $amide[$j], $hamide[$j];
	}
}
close OUP;
