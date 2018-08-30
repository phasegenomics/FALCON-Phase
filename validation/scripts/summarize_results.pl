#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

summarize_results.txt processed_results.txt

Description:

Feed in the table from process_vcf.pl

";


my ($help);
my $opt_success = GetOptions('help'    => \$help,
			      );

die $usage if $help || ! $opt_success;

my $file = shift;
die $usage unless $file;
open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

my $n_unknown = 0;
my $l_unknown = 0;
my $l_count   = 0;
my $total_len = 0;

my %DAT;

LINE: while (<$IN>) {
    chomp;
    $l_count++;
    my @line = split /\t/, $_;
    $total_len += $line[-1];
    if($line[-2] eq "unknown"){
	$n_unknown++;
	$l_unknown += $line[-1];
	next LINE;
    }
    if(! defined $DAT{$line[-2]}{$line[-4]}){
	$DAT{$line[-2]}{$line[-4]}{"len"}           = $line[-1];
	$DAT{$line[-2]}{$line[-4]}{"snp_correct"}   = $line[-5];
	$DAT{$line[-2]}{$line[-4]}{"snp_incorrect"} = $line[-6];
    }
    else{
	$DAT{$line[-2]}{$line[-4]}{"len"}           += $line[-1];
	$DAT{$line[-2]}{$line[-4]}{"snp_correct"}   += $line[-5];
	$DAT{$line[-2]}{$line[-4]}{"snp_incorrect"} += $line[-6];
    }
}

print STDOUT "\n\n          ---  REPORT ---      \n\n";
print STDOUT "n_scored:$l_count total_len:$total_len n_unknown:$n_unknown len_unknown:$l_unknown\n\n";

foreach my $key (keys %DAT){
    foreach my $key2 (keys %{$DAT{$key}}){
	print STDOUT join "\t", ("FP_assignment:$key", "Trio_assignment:$key2", "length:$DAT{$key}{$key2}{len}", $DAT{$key}{$key2}{"snp_correct"}, $DAT{$key}{$key2}{"snp_incorrect"}, "\n" );
    }
}

print STDOUT "\n\n";

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------


