#!/usr/bin/env python

# Sarah B. Kingan
# 20 September 2017
# 
# Pacific Biosciences
# Applications Lab
#
# Convert nucmer.coords files for all haplotigs aligned
# to primary contigs into ncbi placement file
#
###########################################################################################################

# import libraries
import numpy as np
import subprocess
import argparse

###########################################################################################################

desc ='''Make haplotig placement file from individual coords file of all haplotigs to primary'''
parser = argparse.ArgumentParser(description=desc)
parser.add_argument("infile", help="coords file")
args = parser.parse_args()

def gap (j, k):
	if j <= k:
		return 0
	else:
		g = int(hd[j,0]) - int(hd[k,1])
		if g < 0:
			return 0
		else:
			return g


# input file
coords_file = args.infile

# read file	
d=np.genfromtxt(open(coords_file, "rb"), dtype="str")
#coords=[pstart_1 pend hstart_1 hend pAlnL hAlnL PID pL hL pID hID]

# unique set of haplotigs
haplotigs=list(set(list(d[:,10])))
# loop through haplotigs, subset array for each
for h in haplotigs:
	output=['hID','hL', 'hStart_0', 'hEnd', 'ori', 'pID', 'pL', 'pStart_0', 'pEnd', 'match', 'alnL', 'mapQ']

	# get data for a haplotig
	subset = []
	for i in range(0,d.shape[0]):
		if d[i,10] == h:
			subset.append(i)	
	n = len(subset)
	hd = d[subset]

	#initialize results array
	results = [{'start_index':-1, 'end_index':-1, 'score': float("-inf")} for l in range(n)]

	# score chained alignments	
	for i in range(0,n):
		score = 0
		for j in range(i,n):
			score += int(hd[j,4]) * float(hd[j,6]) * 0.01 - gap(j, j-1)
			result_index = j - 1
			if results[result_index]['score'] < score:
				results[result_index]['start_index'] = i
				results[result_index]['end_index'] = j
				results[result_index]['score'] = score

	# get best scoring chained alignment
	max = float("-inf")
	best = {'start_index':-1, 'end_index':-1, 'score': float("-inf")}
	for i in range(0, len(results)):
		if results[i]['score'] > max:
			best = results[i]
			max = best['score']

	# load data
	output[0]=hd[0,10]	#htig name
	output[5]=hd[0,9]	#pcontig name
	output[6]=hd[0,7]	#pcontig length
	output[1]=hd[0,8]	#hcontig length
	hstart=int(hd[best['start_index'],2]) # hstart
	hend=int(hd[best['end_index'],3]) # hend
	pstart=int(hd[best['start_index'],0]) # pstart
	pend=int(hd[best['end_index'],1]) # pend


	# start and end coords
	output[2]=hstart-1 # hstart
	output[3]=hend # hend
	output[7]=pstart-1 # pstart
	output[8]=pend # pend

	# orientation
	output[4]='+' # ori
	hori = '+'
	pori = '+'
	if hstart > hend:
		hori = '-'
		output[2]=hend-1
		output[3]=hstart # hend
	if pstart > pend:
		pori = '-'
		output[7]=pend-1 # pstart
		output[8]=pstart # pend

	if hori != pori:
		output[4]='-'

	# aln details, not real data
	output[10]=output[3]-output[2] # alnL
	output[9]=output[10] # match
	output[11]='60' # mapQ		
# print to stdout
	print('\t'.join(str(o) for o in output))
