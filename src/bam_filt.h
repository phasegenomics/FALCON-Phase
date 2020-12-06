
#ifndef BAM_FILT_H
#define BAM_FILT_H


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"



enum counts {PASS  = 0,
	     SAME_SEQID  = 2,
	     LOW_MQ      = 4,
	     XA_SA       = 8,
	     NM          = 16,
	     SMALLCONTIG = 32,
	     EXCLUDE     = 64,
	     SA_ONLY     = 128,
	     UNMAPPED    = 256,
	     DUPLICATE   = 1024};

struct read_pair {
	bam1_t * read1;
	bam1_t * read2;
} rp;


int filter_bam(char ** argv, int argc);

int parse_length_exclude(bam_hdr_t * header);

int mod_dup_header( bam_hdr_t * dup);


#endif /* BAM_FILT_H */
