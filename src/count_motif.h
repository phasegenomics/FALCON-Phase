#ifndef COUNT_MOTIF
#define COUNT_MOTIF

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <getopt.h>
#include "htslib/faidx.h"
#include "dna.h"



struct sequence{
  char     * name;
  int        slen;
  uint32_t   cutcount;
};

struct sequenceInfo{
    int nseq;
    struct sequence * dat;
    faidx_t * fai;
};



/**
 * [count_motif description]
 * @param  motif      [description]
 * @param  motif_size [description]
 * @param  fai        [description]
 * @param  seqName    [description]
 * @param  seq_idx    [description]
 * @return            [description]
 */
int32_t count_cutsite(char * motif,   faidx_t * fai, const char * seqName, int seq_idx);

/**
 * [count_motif_runner description]
 * @param  argv [description]
 * @param  argc [description]
 * @return      [description]
 */
int count_cutsite_runner(char ** argv, int argc);

int print_sequenceInfo(FILE * stream, struct sequenceInfo * sinfo);

struct sequenceInfo * destroy_sequence_info(struct sequenceInfo * seqInfo);

struct sequenceInfo * load_seq_info(char * fasta, char * motif);

#endif
