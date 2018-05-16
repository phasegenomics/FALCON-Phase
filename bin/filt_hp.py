#!/usr/bin/env python

# import libraries
from argparse import ArgumentParser
import numpy as np


# define command line arguments, program description, and help
desc ='''Remove nested haplotigs.'''
parser = ArgumentParser(description=desc)
parser.add_argument("placement", help="haplotigPlacementFile")
args = parser.parse_args()


# get filenames
inFile = args.placement

# functions
def nested(a,b):
        # a nested inside b
        nested = 0
        astart = int(a[7])
        aend = int(a[8])
        bstart = int(b[7])
        bend = int(b[8])
        if astart < bend and astart > bstart and aend > bstart and aend < bend:
                nested = 1
        return nested

# functions
def dup(a,b):
        # a duplicate of b
        dup = 0
        astart = int(a[7])
        aend = int(a[8])
        bstart = int(b[7])
        bend = int(b[8])
        if astart == bstart and aend == bend:
                dup = 1
        return dup


# open file
hp=np.genfromtxt(open(inFile, "rb"), dtype="str")

# list of primary contigs
primaries=list(set(list(hp[:,5])))
primaries.sort()

# list of booleans indicating filtered or not
filtHaplotigs = [0] * hp.shape[0]


# loop through each primary and process all its haplotigs
for p in primaries:

# list of array indices for primary contig
    subset = []
    for i in range(0,hp.shape[0]):
        if hp[i,5] == p:
            subset.append(i) # list of indices in hp for primary contig
# pairwise comparison of each haplotig
    for i in range(len(subset)):
        for j in range(i):
            if nested(hp[subset[i],:],hp[subset[j],:]):
                filtHaplotigs[subset[i]] = 1
            if nested(hp[subset[j],:],hp[subset[i],:]):
                filtHaplotigs[subset[j]] = 1
            if dup(hp[subset[i],:],hp[subset[j],:]):
                filtHaplotigs[subset[j]] = 1 # remove j
# print new file
for k in range(len(filtHaplotigs)):
    if filtHaplotigs[k] == 0:
        print('\t'.join(str(o) for o in hp[k,:]))
