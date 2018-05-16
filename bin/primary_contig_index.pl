#!/usr/bin/perl -w

####################################################################################################
#
#		Sarah B. Kingan
#		Pacific Biosciences
#		26 January 2018
#
#		Title: primary_contig_index.pl
#
#		Project: phase unzip
#	
#		Input: 	1. AB_pairs.txt
#			2. *.fai of minced unzip genome
#
#		Output: list of minced contigs for each primary contig and A:B pairs
#			
#
####################################################################################################

use strict;
use POSIX qw(ceil);

my $usage = "primary_contig_index.pl AB_pairs.txt mincedUnzip.fa.fai\n";

# AB pairs file
my $AB_pairs_file = shift(@ARGV) or die $usage;

# fai
my $fai_file = shift(@ARGV) or die $usage;

## hash for indexing ##
# key = contig name
# value = index
open (FAI, $fai_file);
my %fai_hash;
my $index = 0;
my @line_array;
while (my $line = <FAI>) {
        chomp$line;
        @line_array = split("\t", $line);
        $fai_hash{$line_array[0]} = $index;
	$index++;
}
$fai_hash{$line_array[0]} = $index;


## hash with list of contigs for each primary ##
# key = primary contig name
# value = string with list of contigs indices
my %mincedSet_hash;
foreach my $key (sort keys %fai_hash) {
#	print $key, "\n";
#	print getPrimary($key), "\n";
#	print getIndex($key), "\n";
	my $primary = getPrimary($key);
	my $i = getIndex($key);
	if ( exists $mincedSet_hash{$primary} ) {
		$mincedSet_hash{$primary} .= ",".$i;
# name		$mincedSet_hash{$primary} .= ",".$key;
	}
	else {
		$mincedSet_hash{$primary} = $i;
# name		$mincedSet_hash{$primary} = $key;
	}
}


## hash of haplotig pairs ##
# key = primary contig name
# value = string with list of pairs names (e.g. 'A1:B1,A2:B2,')
my %AB_hash;
open (AB, $AB_pairs_file);
while (my $line = <AB>) {
	if ($line =~ /[0-9]F/) { # ignore header lines
	        chomp$line;
       		 @line_array = split("\t", $line);
		my $A_name = $line_array[0];
		my $B_name = $line_array[1];
		my $A_index = getIndex($A_name);
		my $B_index = getIndex($B_name);
		my $primary = getPrimary($A_name);
	#	print join(",", ($A_name, $B_name, $A_index, $B_index)), "\n";
		if (exists $AB_hash{$primary} ) {
		        $AB_hash{$primary} .= ",".$A_index.":".$B_index;
	# name	        $AB_hash{$primary} .= ",".$A_name.":".$B_name;
	        }
		else {
			$AB_hash{$primary} .= $A_index.":".$B_index;
	# name		$AB_hash{$primary} .= $A_name.":".$B_name;
		}
	}
}
## get ready to print ##
# unique list of primary contigs
my @tmp = keys %mincedSet_hash;
push ( @tmp, keys %AB_hash );
my @primary_list = uniq(@tmp);


foreach my $contig ( @primary_list ) {
	print $contig, "\t";
	if ( exists $mincedSet_hash{$contig} )	{
		print $mincedSet_hash{$contig}, "\t";
	}
	else {
		print "NA", "\t";
	}
	if ( exists $AB_hash{$contig} )  {
                print $AB_hash{$contig}, "\n";
        }
        else {
                print "NA", "\n";
        }
}




## extract primary contig name from any contig ID
sub getPrimary {
	my ($contig) = @_; 
	my $base = 'error';
	if ( $contig =~ /([0-9]{6}F)/ ) {
		$base = $1;
	}
	return $base;
}

## convert contig name to index
sub getIndex {
        my ($name) = @_;
	my $index = 'error';
        if ( exists $fai_hash{$name} ) {
                $index = $fai_hash{$name};
        }
        return $index;
}

sub uniq {
	my @input = @_;
	my %hash;
	foreach my $a (@input) {
		$hash{$a} = '1';
	}
	my @output = (sort { $a cmp $b } keys %hash);
	return @output;
}
