#!/usr/bin/env python
# vim: ts=4 noet:
# python2.7 and python3 should both work. ~cdunn
####################################################################################################
#
#		Sarah B. Kingan
#		Pacific Biosciences
#		18 April 2018 (python version in September)
#
#		Title: emit_haplotigs.py
#
#		Project: FALCON-phase
#	
#		Input: 	phase.txt, BC.bed, clean_unzip_ref_p_h.fa bedtools_path output_format
#
#		Output: output results in fasta and bed format
#			
#
####################################################################################################

###################
### INPUT FILES ###
###################

# phased.txt
#000000F 000000F_001:0-26421 000000F:10327801-10353990 -nan 0.1852 0.0238 0 6
#000000F 000000F_003:0-35964 000000F:21633239-21669154 -nan 0.0375 0.0375 2 7
#000000F 000000F_002:0-43639 000000F:35315359-35358662 -nan 0.1186 0.0172 1 8
#000001F 000001F_003:0-92432 000001F:5464634-5557226 -nan 0.0725 0.0613 5 9

# BC.bed
#000000F 0	10327801
#000000F 10327801	10353990
#000000F 10353990	21633239
#000000F 21633239	21669154

# clean_unzip_asm_p_h.fa
# path_to_bedtools
# output_format

###########################
### DEFINE OUTPUT FILES ###
###########################

bed0 = 'b0.bed';
bed1 = 'b1.bed';
fasta0 = 'f0.fa';
fasta1 = 'f1.fa';
tab0 = 'tmp_phase0.txt';
tab1 = 'tmp_phase1.txt';

import collections, logging, os, sys
LOG = logging.getLogger()

PhasePair = collections.namedtuple('PhasePair', 'phase pair')

def read_phase_file(stream):
	###############################
	### PHASE HASH FOR AB PAIRS ###
	# ID
	# 	phase = [0,1]
	#	pair = ID of other member of pair
	###############################
	phase_hash = collections.defaultdict(PhasePair)
	for line in stream:
		line = line.rstrip()
		line_array = line.split()
		if line_array[3] == '-nan':
			LOG.warning("error in phasing, did you run enough iterations ('iter' >= 10e6)? {} and {} 'phased.txt' contains '-nan'".format(
			line_array[1], line_array[2]))
		phase_hash[line_array[1]] = PhasePair(phase=0, pair=line_array[2])
		phase_hash[line_array[2]] = PhasePair(phase=1, pair=line_array[1])
	return phase_hash

def main(prog, phase_file, BC_bed_file, unzip_asm_file, bedtools_path, output_format):
	#my $usage = "emit_haplotigs.pl phased.txt BC.bed clean_unzip_asm_p_h.fa path_to_bedtools output_format\n";
	#if ($output_format eq 'unzip') {
	#	$bed0 = 'cns_p_ctg_phased.bed';
	#	$bed1 = 'cns_h_ctg_phased.bed';
	#	$fasta0 = 'cns_p_ctg_phased.fa';
	#	$fasta1 = 'cns_h_ctg_phased.fa';
	#}
	#elsif ($output_format eq 'pseudohap') {
	#	# do nothing
	#}
	logging.basicConfig()
	assert output_format in ('unzip', 'pseudohap'), "Please specify output format: 'pseudohap' or 'unzip'"
	with open(phase_file) as stream:
		global phase_hash
		phase_hash = read_phase_file(stream)
	##########################
	#### PHASED BED FILES ####
	##########################
	makeBED(BC_bed_file, output_format);
	######################
	#### RUN BEDTOOLS ####
	######################
	# use bedtools to make tab delimited seq files
	#unless (-s $phase0_tab_file)
	cmd = "{} getfasta -tab -name -fi {} -bed {} > {}".format(
			bedtools_path, unzip_asm_file, bed0, tab0)
	syscall(cmd);
	#unless (-s $phase1_tab_file)
	cmd = "{} getfasta -tab -name -fi {} -bed {} > {}".format(
			bedtools_path, unzip_asm_file, bed1, tab1)
	syscall(cmd);
	########################
	#### GENERATE FASTA ####
	########################
	prefix0 = '';
	prefix1 = '';
	if (output_format == 'pseudohap'):
		prefix0 = '_0';
		prefix1 = '_1';
	tab2fa(tab0, prefix0, fasta0);
	tab2fa(tab1, prefix1, fasta1);
	### TODO: add option for no suffix for unzip fmt

def syscall(cmd, unchecked=False):
	rc = os.system(cmd)
	if not unchecked:
		if rc:
			msg = '{} <- {!r}'.format(rc, cmd)
			raise Exception(msg)
	return rc

def get_hID(IDa, IDb):
	"""
	>>> get_hID('ab', 'cd')
	''
	>>> get_hID('a_b', 'cd')
	'a_b'
	>>> get_hID('ab', 'c_d')
	'c_d'
	"""
	if ('_' in IDa):
		hID = IDa
	elif ('_' in IDb):
		hID = IDb;
	else:
		hID = '';
	return hID;

"""
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
"""

def switch2(zero, one):
	"""
	>>> switch2('', '_')
	False
	>>> switch2('_', '')
	True
	"""
	zero_count = zero.count('_')
	one_count = one.count('_')
	return (zero_count > one_count)

def _tab2fa(TAB, suffix, FA):
	ID = 'first';
	seq = list()
	for line in TAB:
		line = line.rstrip()
		line_array = line.split('\t')
		if (line_array[0] == ID):
			seq.append(line_array[1])
		else:
			if ID != 'first':
				FA.write('>{}{}\n'.format(ID, suffix))
				FA.write('{}\n'.format(''.join(seq)))
			ID = line_array[0]
			seq = [line_array[1]]
	FA.write('>{}{}\n'.format(ID, suffix))
	FA.write('{}\n'.format(''.join(seq)))

def tab2fa(in_fn, suffix, out_fn):
	with open(in_fn) as TAB, open(out_fn, 'w') as FA:
		_tab2fa(TAB, suffix, FA)

def _makeBED_unzip(BC, ZERO, ONE):
	zero_bed = list()
	one_bed = list()
	oldPID = 'first'
	for line in BC:
		line = line.rstrip()
		line_array = line.split('\t')
		contig_ID_BC = '{}:{}-{}'.format(*line_array)
		primaryID = line_array[0]
		tup = phase_hash.get(contig_ID_BC, None)
		if tup:
			if tup.phase == 0:
				zero_entry = '{}\t{}'.format(line, primaryID)
				contig_ID_A = tup.pair.split(':-')
				haplotigID = contig_ID_A[0]
				contig_ID_A.append(haplotigID)
				one_entry = '\n'.join(contig_ID_A)
			else:
				one_entry = '{}\t{}'.format(line, primaryID)
				contig_ID_A = tup.pair.split(':-')
				haplotigID = contig_ID_A[0]
				contig_ID_A.append(haplotigID)
				zero_entry = '\n'.join(contig_ID_A)
		else: # collapsed haplotype
			zero_entry = '{}\t{}.collapsed'.format(line, primaryID)
			one_entry = '{}\t{}.collapsed'.format(line, primaryID)
		if (primaryID != oldPID and oldPID != 'first'): # new primary
				if (switch2(''.join(zero_bed), ''.join(one_bed))): # 0 phase most resembles original primary contig
					zero_bed, one_bed = one_bed, zero_bed
				for i in range(len(zero_bed)):
					if ('collapsed' in zero_bed[i]): # clean it
						zero_bed[i] = zero_bed[i].replace('collapsed', '', 1)
					else:
						if ('_' in zero_bed[i]): # swap names so haplotigs IDs are in h file
							array0 = zero_bed[i].split('\t')
							array1 = one_bed[i].split('\t')
							array0[3], array1[3] = array1[3], array0[3]
							zero_bed[i] = '\t'.join(array0)
							one_bed[i] = '\t'.join(array1)
						ONE.write('{}\n'.format(one_bed[i]))
					ZERO.write('{}\n'.format(zero_bed[i]))
				zero_bed = list()
				one_bed = list()
		zero_bed.append(zero_entry)
		one_bed.append(one_entry)
		oldPID = primaryID
	if (switch2(''.join(zero_bed), ''.join(one_bed))): # 0 phase most resembles original primary contig
		zero_bed, one_bed = one_bed, zero_bed
	for i in range(len(zero_bed)):
		if 'collapsed' in zero_bed[i]:
			zero_bed[i] = zero_bed[i].replace('collapsed', '', 1)
		else:
			if '_' in zero_bed[i]: # swap names so haplotigs IDs are in h file
				array0 = zero_bed[i].split('\t')
				array1 = one_bed[i].split('\t')
				array0[3], array1[3] = array1[3], array0[3]
				zero_bed[i] = '\t'.join(array0)
				one_bed[i] = '\t'.join(array1)
			ONE.write('{}\n'.format(one_bed[i]))
		ZERO.write('{}\n'.format(zero_bed[i]))

def _makeBED_pseudohap(BC, ZERO, ONE):
	for line in BC:
		line = line.rstrip()
		line_array = line.split('\t')
		contig_ID = '{}:{}-{}'.format(*line_array)
		tup = phase_hash.get(contig_ID, None)
		if tup:
			if tup.phase == 0:
				zero_bed = line
				one_bed = tup.pair.replace(':', '\t').replace('-', '\t')
			else:
				assert tup.phase == 1, tup
				one_bed = line
				zero_bed = tup.pair.replace(':', '\t').replace('-', '\t')
		else: # collapsed
			zero_bed = line
			one_bed = line
		ZERO.write('{}\t{}\n'.format(zero_bed, line_array[0]))
		ONE.write('{}\t{}\n'.format(one_bed, line_array[0]))

def makeBED(BC_fn, fmt):
	with open(BC_fn) as BC, open(bed0, 'w') as ZERO, open(bed1, 'w') as ONE:
		func = {
			'unzip': _makeBED_unzip,
			'pseudohap': _makeBED_pseudohap,
		}[fmt]
		func(BC, ZERO, ONE)

def test():
	import doctest
	doctest.testmod()

if __name__ == "__main__":
	if '--test' in sys.argv:
		sys.exit(test())
	main(*sys.argv)
