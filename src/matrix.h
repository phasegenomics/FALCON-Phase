
#ifndef MATRIX_H
#define MATRIX_H

#ifndef MAGIC_HEAD
#define MAGIC_HEAD 1984
#endif

#ifndef MAGIC_TAIL
#define MAGIC_TAIL 1973
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
allows the matrix to be an arbitrary type
 */
typedef uint32_t datum;
typedef  double fwf;

// the matrix struct
struct matrix {
    datum n1;
    datum n2;
    fwf ** dat;
};

/**
 * allocates the memory for the matrix and sets all values equal to zero
 * @param  n1 - the number for x dimension
 * @param  n2 - the number for y dimension
 * @return   - a matrix pointer
 */
struct matrix * init_matrix(const datum n1, const datum n2);

/**
 * frees the memory for the matrix;
 * @param matrix - the matrix to destroy
 */
void destroy_matrix(struct matrix * m);

/**
 * Counts into the matrix, checks bounds
 * @param  m - matrix pointer
 * @param  x - index 1
 * @param  y - index 2
 * @return   0 if bounds are bad ; 1 if good
 */
int add_link(struct matrix * m, datum x, datum y);

/**
 * Prints the full matrix
 * @param m  - matrix pointer
 * @param in - 1 if sparse ; full
 */
int print_matrix(struct matrix * m, const char * fn_out);

/**
 * serialize the 2D matrix
 * @param m   - matrix pointer
 * @param fn  - filename to write to
 * @return < 1 if failure
 */
int freeze_matrix(struct matrix * m, const char * fn_out);

/**
 * reads serialized matrix from disk to 2D matrix
  * @param  fn_in - file to read in
 * @return NULL if failed to thaw matrix
 */
struct matrix * thaw_matrix(const char * fn_in);

/**
 * This function sums the matrix
 * @param  mat - the matrix
 * @return  datum the sum of the matrix
 */
datum matrix_sum(struct matrix * mat);

int print_matrix_bounded(struct matrix * m, datum i_max, datum j_max);

int print_matrix_bounded_square(struct matrix * m, datum i_max, datum j_max);

#endif /* MATRIX_H */
