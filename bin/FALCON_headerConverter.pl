#!/usr/bin/perl -w

####################################################################################################
#
#		Sarah B. Kingan
#		Pacific Biosciences
#		15 August 2018
#
#		Title: FALCON_headerConverter.pl
#
#		Project: FALCON-Phase
#	
#		Input: 	FASTA files for FALCON (post March 2018)
#
#		Output: FASTA file with "old style" FALCON headers to STDOUT
#				
#			
#
####################################################################################################

use strict;

my $usage = "FALCON_headerConverter.pl cns_X_ctg.fasta\n";

my $in = shift(@ARGV) or die $usage;
#000877Fp01_00001

my %p_hash;

open (IN, $in);
while (my $line = <IN>) {
	if ($line !~ /^>/) {
		print $line;
	}
	elsif ($line =~ /^(>[0-9]{6}F)(p*[0-9]*)_([0-9]+)(\|arrow)/) { # haplotig
		$p_hash{$2} = 1;
		print $1, "_";
		print sprintf("%03d",$3);
		print $4, "\n";
	}
	elsif ($line =~ /^(>[0-9]{6}F)(p[0-9]+)(\|arrow)/) { # primary contig
		$p_hash{$2} = 1;
		print $1;
		print $3, "\n";
	}
}

if (scalar(keys%p_hash) > 1) {
	print STDERR "ERROR: Conversion will result in loss of contigs in FALCON-Phase. You have more than one contig with same base name [0-9]{6} but different p[0-9]+ suffix. Identify these and rename them manually.\n";
}
