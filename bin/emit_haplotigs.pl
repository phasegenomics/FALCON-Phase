#!/usr/bin/perl -w

####################################################################################################
#
#		Sarah B. Kingan
#		Pacific Biosciences
#		18 April 2018
#
#		Title: emit_haplotigs.pl
#
#		Project: FALCON-phase
#	
#		Input: 	phase.txt, BC.bed, clean_unzip_ref_p_h.fa
#
#		Output: phased haplotigs in fasta, plus bed of minced contigs in each emitted haplotig
#			
#
####################################################################################################

use strict;
use warnings;

my $usage = "emit_haplotigs.pl phased.txt BC.bed clean_unzip_asm_p_h.fa pwd/bedtools\n";


###################
### INPUT FILES ###
###################

my $phase_file = shift(@ARGV) or die $usage;
#000000F 000000F_001:0-26421 000000F:10327801-10353990 -nan 0.1852 0.0238 0 6
#000000F 000000F_003:0-35964 000000F:21633239-21669154 -nan 0.0375 0.0375 2 7
#000000F 000000F_002:0-43639 000000F:35315359-35358662 -nan 0.1186 0.0172 1 8
#000001F 000001F_003:0-92432 000001F:5464634-5557226 -nan 0.0725 0.0613 5 9

my $BC_bed_file = shift(@ARGV) or die $usage;
#000000F 0       10327801
#000000F 10327801        10353990
#000000F 10353990        21633239
#000000F 21633239        21669154

my $unzip_asm_file = shift(@ARGV) or die $usage;

my $bedtools_path  = shift(@ARGV) or die $usage;

####################
### OUTPUT FILES ###
####################
my $phase0_bed_file = 'phase0.bed';
my $phase1_bed_file = 'phase1.bed';
my $phase0_tab_file = 'tmp_phase0.txt';
my $phase1_tab_file = 'tmp_phase1.txt';
my $fasta_file = 'diploid_phased.fa';


###############################
### MINCED CONTIGS AB PAIRS ###
# ID
# 	phase = [0,1]
#	pair = ID
###############################
# make phase hash
my %phase_hash;
open (PH, $phase_file) or die "Could not open file '$phase_file'";
while (my $line = <PH>) {
	chomp $line;
	my @line_array = split(" ", $line);
	if ($line_array[3] eq "-nan") {
		print STDERR "Warning: error in phasing ", $line_array[1], " and ", $line_array[2], " 'phased.txt' contains '-nan'\n";
	}
	$phase_hash{$line_array[1]}{'phase'} = 0;
	$phase_hash{$line_array[2]}{'phase'} = 1;
	$phase_hash{$line_array[1]}{'pair'} = $line_array[2];
	$phase_hash{$line_array[2]}{'pair'} = $line_array[1];
}
close $phase_file;

##########################
#### PHASED BED FILES ####
##########################
# create phase 0 and phase 1 bed files
open (BC, $BC_bed_file) or die "Could not open file '$BC_bed_file'";
open (ZERO, '>', $phase0_bed_file) or die "Could not open file '$phase0_bed_file'";
open (ONE, '>', $phase1_bed_file) or die "Could not open file '$phase1_bed_file'";
while (my $line = <BC>) {
	chomp $line;
	my @line_array = split("\t", $line);
	my $contig_ID = $line_array[0].":".$line_array[1]."-".$line_array[2];
	my $zero_bed;
	my $one_bed;
#	print $contig_ID, "\n";
	if (exists $phase_hash{$contig_ID}) {
		if ($phase_hash{$contig_ID}{'phase'} == 0) {
			$zero_bed = $line;
			$one_bed = $phase_hash{$contig_ID}{'pair'};
			$one_bed =~ s/:/\t/;
			$one_bed =~ s/-/\t/;
		}
		elsif ($phase_hash{$contig_ID}{'phase'} == 1) {
			$one_bed = $line;
			$zero_bed = $phase_hash{$contig_ID}{'pair'};
			$zero_bed =~ s/:/\t/;
			$zero_bed =~ s/-/\t/;
		}
	}		
	else { # collapsed
		$zero_bed = $line;
		$one_bed = $line;
	}
	print ZERO $zero_bed, "\t", $line_array[0], "\n";
	print ONE $one_bed, "\t", $line_array[0], "\n";
}


########################
#### GENERATE FASTA ####
########################
# use bedtools to make tab delimited seq files
my $cmd;
#unless (-s $phase0_tab_file) {
	$cmd = "$bedtools_path getfasta -tab -name -fi $unzip_asm_file -bed $phase0_bed_file > $phase0_tab_file";
	system($cmd);
#}
#unless (-s $phase1_tab_file) {
	$cmd = "$bedtools_path getfasta -tab -name -fi $unzip_asm_file -bed $phase1_bed_file > $phase1_tab_file";
	system($cmd);
#}
# process into concatenated fasta
my $seq;
my $pcontig;
open (P0T, $phase0_tab_file) or die "Could not open file '$phase0_tab_file'";
#open (FA, '>', $fasta_file) or die "Could not open file '$fasta_file'";
$pcontig = 'first';
$seq = '';
while (my $line = <P0T>) {
	chomp $line;
	my @line_array = split("\t", $line);
	if ($line_array[0] eq $pcontig) {
		$seq .= $line_array[1];
	}
	else {
		if ($pcontig eq 'first') {
			# don't print
		}
		else {
			print ">", $pcontig, "_0", "\n";
			print $seq, "\n";
		}
		$pcontig = $line_array[0];
		$seq = $line_array[1];
	}
}
print ">", $pcontig, "_0", "\n";
print $seq, "\n";
close $phase0_tab_file;

open (P1T, $phase1_tab_file) or die "Could not open file '$phase1_tab_file'";
$pcontig = 'first';   
$seq = '';
while (my $line = <P1T>) {
        chomp $line;
        my @line_array = split("\t", $line);
        if ($line_array[0] eq $pcontig) {
                $seq .= $line_array[1];
        }
        else {
                if ($pcontig eq 'first') {
                        # don't print
                }
                else {
                        print ">", $pcontig, "_1", "\n";
                        print $seq, "\n";
                }
                $pcontig = $line_array[0];
                $seq = $line_array[1];
        }
}
print ">", $pcontig, "_1", "\n";
print $seq, "\n";
close $phase1_tab_file;
#close $fasta_file;




exit;
