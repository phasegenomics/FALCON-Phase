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
#		Input: 	phase.txt, BC.bed, clean_unzip_ref_p_h.fa bedtools_path output_format
#
#		Output: output results in fasta and bed format
#			
#
####################################################################################################

use strict;
use warnings;

my $usage = "emit_haplotigs.pl phased.txt BC.bed clean_unzip_asm_p_h.fa path_to_bedtools output_format\n";


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

my $output_format = shift(@ARGV) or die $usage;

###########################
### DEFINE OUTPUT FILES ###
###########################

my $bed0 = 'b0.bed';
my $bed1 = 'b1.bed';
my $fasta0 = 'f0.fa';
my $fasta1 = 'f1.fa';
my $tab0 = 'tmp_phase0.txt';
my $tab1 = 'tmp_phase1.txt';
#if ($output_format eq 'unzip') {
#	$bed0 = 'cns_p_ctg_phased.bed';
#	$bed1 = 'cns_h_ctg_phased.bed';
#	$fasta0 = 'cns_p_ctg_phased.fa';
#	$fasta1 = 'cns_h_ctg_phased.fa';
#}
#elsif ($output_format eq 'pseudohap') {
#	# do nothing
#}
unless (($output_format eq 'unzip') || ($output_format eq 'pseudohap')) {
	print STDERR "Please specify output format: 'pseudohap' or 'unzip'\n";
	die;
}


###############################
### PHASE HASH FOR AB PAIRS ###
# ID
# 	phase = [0,1]
#	pair = ID of other member of pair
###############################
my %phase_hash;
open (PH, $phase_file) or die "Could not open file '$phase_file'";
while (my $line = <PH>) {
	chomp $line;
	my @line_array = split(" ", $line);
	if ($line_array[3] eq "-nan") {
		print STDERR "Warning: error in phasing, did you run enough iterations ('iter' >= 10e6)? ", $line_array[1], " and ", $line_array[2], " 'phased.txt' contains '-nan'\n";
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
makeBED($BC_bed_file, $output_format);


######################
#### RUN BEDTOOLS ####
######################
# use bedtools to make tab delimited seq files
my $cmd;
#unless (-s $phase0_tab_file) {
	$cmd = "$bedtools_path getfasta -tab -name -fi $unzip_asm_file -bed $bed0 > $tab0";
	system($cmd);
#}
#unless (-s $phase1_tab_file) {
	$cmd = "$bedtools_path getfasta -tab -name -fi $unzip_asm_file -bed $bed1 > $tab1";
	system($cmd);
#}


########################
#### GENERATE FASTA ####
########################
my $prefix0 = '';
my $prefix1 = '';
if ($output_format eq 'pseudohap') {
	$prefix0 = '_0';
	$prefix1 = '_1';
}
tab2fa($tab0, $prefix0, $fasta0);
tab2fa($tab1, $prefix1, $fasta1);


### add option for no suffix for unzip fmt

exit;

sub switch1 {
	my (@zero_bed_array) = @_;
	my $count = 0;
	my $half = scalar(@zero_bed_array)/ 2;
	my $switch = 0;
	for (my $i = 0; $i < scalar(@zero_bed_array); $i++) {
		if ($zero_bed_array[$i] =~ /_/) {
			$count++;
		}
	}
	if ($count > $half) {
		$switch = 1;
	}
	return $switch;
}

sub get_hID {
	my ($IDa, $IDb) = @_;
	my $hID;
	if ($IDa =~ '_') {
		$hID = $IDa;
	}
	elsif ($IDb =~ '_') {
		$hID = $IDb;
	}
	return $hID;
}

sub switch2 {
	my ($zero, $one) = (@_);
	my $switch = 0;
	my $zero_count = ($zero =~ tr/_//);
	my $one_count = ($one =~ tr/_//);
	if ($zero_count > $one_count) {
		$switch = 1;
	}
	return $switch;
}

sub tab2fa {
	my ($in, $suffix, $out) = (@_);
	open(TAB, $in) or die "Could not open file '$in'";
	open(FA, ">", $out) or die "Could not open file '$out'";
	my $ID = 'first';
	my $seq = '';
	while (my $line = <TAB>) {
		chomp $line;
		my @line_array = split("\t", $line);
		if ($line_array[0] eq $ID) {
			$seq .= $line_array[1];
		}
		else {
			unless ($ID eq 'first') {
				print FA ">", "$ID", $suffix, "\n";
				print FA $seq, "\n";
			}
			$ID = $line_array[0];
			$seq = $line_array[1];
		}
	}
	print FA ">", "$ID", $suffix, "\n";
	print FA $seq, "\n";
	close $in;
	close $out;
}


sub makeBED {
	my ($BC, $fmt) = (@_);
	open (BC, $BC) or die "Could not open file '$BC'";
	open (ZERO, '>', $bed0) or die "Could not open file '$bed0'";
	open (ONE, '>', $bed1) or die "Could not open file '$bed1'";

	if ($fmt eq 'unzip') {
		my @zero_bed;
		my @one_bed;
		my $contig_ID_BC;
		my @contig_ID_A;
		my $oldPID = 'first';
		my $primaryID;
		my $haplotigID;
		my $zero_entry;
		my $one_entry;
		while (my $line = <BC>) {
			chomp $line;
			my @line_array = split("\t", $line);
			$contig_ID_BC = $line_array[0].":".$line_array[1]."-".$line_array[2];
			$primaryID = $line_array[0];
			if (exists $phase_hash{$contig_ID_BC}) { # haplotig block
				if ($phase_hash{$contig_ID_BC}{'phase'} == 0) { # primary phase 0
					$zero_entry = $line."\t".$primaryID;
					@contig_ID_A = split(/:|-/, $phase_hash{$contig_ID_BC}{'pair'});
					$haplotigID = $contig_ID_A[0];
					$one_entry = join("\t", (@contig_ID_A, $haplotigID));
				}
				elsif ($phase_hash{$contig_ID_BC}{'phase'} == 1) {
					$one_entry = $line."\t".$primaryID;
					@contig_ID_A = split(/:|-/, $phase_hash{$contig_ID_BC}{'pair'});
					$haplotigID = $contig_ID_A[0];
					$zero_entry =  join("\t", (@contig_ID_A, $haplotigID));
				}
			}
			else { # collapsed haplotype
				$zero_entry = $line."\t".$primaryID."collapsed";
				$one_entry = $line."\t".$primaryID."collapsed";
			}
			if ($primaryID ne $oldPID) { # new primary
				unless ($oldPID eq 'first') {
					if (switch2(join("",@zero_bed),join("",@one_bed)) == 1) { # 0 phase most resembles original primary contig
						my @tmp = @zero_bed;
						@zero_bed = @one_bed;
						@one_bed = @tmp;
					}
					for (my $i = 0; $i < scalar(@zero_bed); $i++) {
						if ($zero_bed[$i] =~ 'collapsed') { # clean it
							$zero_bed[$i] =~ s/collapsed//;
						}
						else {
							if ($zero_bed[$i] =~ /_/) { # swap names so haplotigs IDs are in h file
								my @array0 = split("\t", $zero_bed[$i]);
								my @array1 = split("\t", $one_bed[$i]);
								my $tmp = $array0[3];
								$array0[3] = $array1[3];
								$array1[3] = $tmp;
								$zero_bed[$i] = join("\t",@array0);
								$one_bed[$i] = join("\t",@array1);
							}
							print ONE $one_bed[$i], "\n";
						}
						print ZERO $zero_bed[$i], "\n";
					}
					@zero_bed = ();
					@one_bed = ();
				}
			}
			push(@zero_bed, $zero_entry);
			push(@one_bed, $one_entry);
			$oldPID = $primaryID;
		}
		if (switch2(join("",@zero_bed),join("",@one_bed)) == 1) { # 0 phase most resembles original primary contig
			my @tmp = @zero_bed;
			@zero_bed = @one_bed;
			@one_bed = @tmp;
		}
		for (my $i = 0; $i < scalar(@zero_bed); $i++) {
			if ($zero_bed[$i] =~ 'collapsed') {
				$zero_bed[$i] =~ s/collapsed//;
			}
			else {
				if ($zero_bed[$i] =~ /_/) { # swap names so haplotigs IDs are in h file
					my @array0 = split("\t", $zero_bed[$i]);
					my @array1 = split("\t", $one_bed[$i]);
					my $tmp = $array0[3];
					$array0[3] = $array1[3];
					$array1[3] = $tmp;
					$zero_bed[$i] = join("\t",@array0);
					$one_bed[$i] = join("\t",@array1);
				}
				print ONE $one_bed[$i], "\n";
			}
			print ZERO $zero_bed[$i], "\n";
		}
	}



	elsif ($fmt eq 'pseudohap') {
		while (my $line = <BC>) {
			chomp $line;
			my @line_array = split("\t", $line);
			my $contig_ID = $line_array[0].":".$line_array[1]."-".$line_array[2];
			my $zero_bed;
			my $one_bed;
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
	}
}
