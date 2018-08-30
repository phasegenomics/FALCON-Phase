#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
my $usage = "

Synopsis:

process_vcf.pl -p 9 -m 10 -d 2,200 your.vcf FP.results.txt

Description:

 -p : VCF column for paternal genotype, zero based, >8 (info column).
 -m : VCF column for maternal genotype, zero based, >8 (info column).
 -d : High and low depth filters.

 VCF cannot be compressed.


";


my ($help);
my $pat;
my $mat;
my $depth;
my $ld;
my $hd;
my $opt_success = GetOptions('help'    => \$help,
			     'p=i' => \$pat,
			     'm=i' => \$mat,
			     'd=s' => \$depth,
    );

die $usage if $help || ! $opt_success;

my $file      = shift;
my $to_score  = shift;

die $usage unless $file && $to_score && $pat && $mat && $depth;

($ld, $hd) = split ",", $depth;

open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";

my %DATA;
my %MDAT;
my %OK;

$OK{"0|0"} = 1;
$OK{"0/0"} = 1;
$OK{"1|1"} = 1;
$OK{"1/1"} = 1;

my $count = 0;

LINE: while (<$IN>) {
    
    chomp;
    next LINE if $_ =~ /^\#/;
    my @line = split /\t/, $_;

    my @mat_dat = split ":", $line[$mat];
    my @pat_dat = split ":", $line[$pat];

    next LINE if(!defined $OK{$mat_dat[0]} || !defined $OK{$pat_dat[0]});
    next LINE if($mat_dat[0] eq $pat_dat[0]);

    if(! defined $DATA{$line[0]}{$mat_dat[0]}){
	$DATA{$line[0]}{$mat_dat[0]} = 1;
    }
    else{
	$DATA{$line[0]}{$mat_dat[0]}++;
    }
    $count ++;
}

print STDERR "SNPs used: $count\n";

close $IN;

foreach my $k (keys %DATA){
    
    my $assignment = "M";

    my $hom_r = 0;
    my $hom_a = 0;

    if(defined $DATA{$k}{"0/0"}){
	$hom_r = $DATA{$k}{"0/0"};
    }
    if(defined $DATA{$k}{"1/1"}){
	$hom_a = $DATA{$k}{"1/1"};
    }
    if($hom_a > $hom_r){
	$assignment = "P";
	my $tmp = $hom_a;
	$hom_a = $hom_r;
	$hom_r = $tmp;
    }

    my $freq = sprintf("%.3f", $hom_r/($hom_r+$hom_a));

    $MDAT{$k} = join "\t", ($k, $hom_r, $hom_a, $hom_r + $hom_a, $assignment, $freq);

}

close $IN;


open (my $INB, '<', $to_score) or die "Can't open $to_score for reading\n$!\n";

LINE: while (<$INB>) {
    chomp;
    my @line = split /\s+/, $_;
    
    my $ans = "unknown";
    my $hap = "unknown";
    my $n;
    my $s;
    my $e;

    ($n, $s, $e) = split /:|-/, $line[1];
    my $l1 = $e - $s;
    ($n, $s, $e) = split /:|-/, $line[1];
    my $l2 = $e - $s;

    if(defined $MDAT{$line[1]}){
	$ans = $MDAT{$line[1]};
	$hap = "hap0";
    }
    if(defined $MDAT{$line[2]}){
	$ans = $MDAT{$line[2]};
	$hap = "hap1";
    }
    
    print "$_\t$ans\t$hap\t$l1\n";

}

close $INB;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

sub {

}

