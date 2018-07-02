#ifndef RUN_PHASING_H
#define RUN_PHASING_H


#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <float.h>

#include "matrix.h"
#include "count_motif.h"
#include "string_parser.h"

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

struct ctg_overlap{
  int first      ;
  int second     ;
  int label_one  ;
  int label_two  ;
  int label_count_one;
  int label_count_two;
};



struct line_data{
  int    n_overlap              ;
  int    n_ctg                  ;
  int    has_overlap            ;
  char   * group                ;
  int    * indices              ;
  struct ctg_overlap * overlaps ;
  struct line_data * before     ;
  struct line_data * after      ;
};

struct phasing_options{
  int       verbose;
  double    nsweeps;
  double     burnin;
  double    damping;
  char  *     fasta;
  char  *    binmat;
  char  *     motif;
  char  *    prefix;
  char  *     index;
}global_opts_phasing;


/**
 * Runs the clustering code
 * @param  argv - command line options
 * @param  argc - command line index
 * @return       > 0 if problem
 */
int run_phasing(char ** argv, int argc );

/**
 * Print the index and value of the vector
 * @param v - the vector
 * @param l - 0-l
 */
void printv(double * v, datum l);

int normalize_matrix_phasing(struct matrix * mat, struct sequenceInfo * seq);

/**
 * Parses the command line options into the global_opts
 * @param  argc [description]
 * @param  argv [description]
 * @return      [description]
 */
int parse_command_line_phasing(char ** argv, int argc);



#endif /* RUN_IRM_H */
