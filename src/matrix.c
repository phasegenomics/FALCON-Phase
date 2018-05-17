/*
Author: Zev Kronenberg
Contact :zev@phasegenomics.com
Date: May 17th 2018

The Clear BSD + Attribution License

Copyright (c) 2018, Pacific Biosciences of California, Inc. and Phase Genomics, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are 
permitted (subject to the limitations in the disclaimer below) provided that the 
following conditions are met:

1.Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2.Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation 
and/or other materials provided with the distribution.

3.All advertising materials mentioning features or use of this software 
must display the following acknowledgement:
  This <product/service> includes <software/the use of software> developed 
  by Pacific Biosciences of California, Inc. and Phase Genomics, Inc.

  4.Distributions of data generated through the use of this software as 
part of a service provided by a for-profit organization must be accompanied 
by the above copyright notice, this list of conditions, the following 
acknowledgement and the disclaimer below:
  This data was generated using software developed by Pacific Biosciences of 
  California, Inc. and Phase Genomics, Inc.

  5.Neither the names of the copyright holders nor the names of their 
contributors may be used to endorse or promote products derived from this 
software without specific prior written permission.

NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY’S PATENT RIGHTS ARE GRANTED BY 
THIS LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS 
OR THEIR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



#include "matrix.h"

datum matrix_sum(struct matrix * mat){

        datum x   = 0;
        datum y   = 0;
        datum sum = 0;

//  fprintf(stderr, "N: %i by N: %i \n", mat->n1, mat->n2 );

        for(; x < mat->n1; x++ ) {
                y = 0;
                for(; y < mat->n2; y++) {
                        sum += mat->dat[x][y];
                }
        }

        return sum;
}

struct matrix * init_matrix(const datum n1, const datum n2)
{
        struct matrix * link_data = malloc(sizeof(struct matrix));

        link_data->n1 = n1;
        link_data->n2 = n2;

        link_data->dat = malloc(sizeof(fwf*)*n1);

        datum i = 0;
        for(; i < n1; i++ ) {
                link_data->dat[i] = malloc(sizeof(fwf)*n2);
                memset( link_data->dat[i], 0, sizeof(fwf)*n2);
        }
        return link_data;
}

void destroy_matrix(struct matrix * m)
{
        datum i = m->n1;
        for(; i < m->n1; i++) {
                free(m->dat[i]);
        }
        free(m);
}

int add_link(struct matrix * m, datum x, datum y)
{
        if (x >= m->n1) return 0;
        if (y >= m->n2) return 0;

        //  fprintf(stderr, "b%i %i %i - ", x, y,  matrix->dat[x][y]);
        m->dat[x][y] += 1;
        //  fprintf(stderr, "b%i %i %i\n ", x, y,  matrix->dat[x][y]);

        return 1;
}

int print_matrix(struct matrix * m, const char * fn_in)
{

        FILE * fh;
        fh = fopen(fn_in, "wb");
        if(fh == NULL) return 1;

        fprintf(stderr, "INFO: matrix_size: %i by %i\n", m->n1, m->n2);
        fprintf(fh, "x\ty\tz\n");

        datum ol = 0;

        for(; ol < m->n1; ol++) {
                datum il = 0;
                for(; il < m->n2; il++) {
                        fprintf(fh, "%i\t%i\t%f\n", ol, il, m->dat[ol][il]);
                }
        }

        fclose(fh);
        return 0;
}


int print_matrix_bounded(struct matrix * m, datum i_max, datum j_max)
{

        fprintf(stderr, "INFO: matrix_size: %i by %i\n", m->n1, m->n2);
        fprintf(stdout, "x\ty\tz\n");

        datum ol = 0;
        datum il = 0;

        for(; ol < i_max; ol++) {
                il = 0;
                for(; il < j_max; il++) {
                        fprintf(stdout, "%i\t%i\t%f\n", ol, il, m->dat[ol][il]);
                }
        }
        return 0;
}

int print_matrix_bounded_square(struct matrix * m, datum i_max, datum j_max)
{

        fprintf(stderr, "INFO: matrix_size: %i by %i\n", m->n1, m->n2);

        datum ol = 0;
        datum il = 0;

        for(; ol < i_max; ol++) {
                il = 0;
                for(; il < j_max; il++) {
                        fprintf(stdout, "%f\t", m->dat[ol][il]);
                }
                fprintf(stdout, "\n");
        }
        return 0;
}


int freeze_matrix(struct matrix * m, const char * fn_out)
{
        FILE * fh;
        fh = fopen(fn_out, "wb");
        if(fh == NULL) return 1;

        uint64_t mh = MAGIC_HEAD;
        uint64_t mt = MAGIC_TAIL;

        fwrite(&mh, sizeof(uint64_t), 1, fh);

        fwrite(&m->n1, sizeof(datum), 1, fh);
        fwrite(&m->n2, sizeof(datum), 1, fh);

        datum index = 0;
        for(; index < m->n1; index += 1)
        {
                fwrite(m->dat[index], sizeof(fwf), m->n2, fh);
        }

        fwrite(&mt, sizeof(uint64_t), 1, fh);

        return 0;
}

struct matrix * thaw_matrix(const char * fn_in ){

        uint64_t magicFront = 0;
        uint64_t magicBack  = 0;

        FILE * fh;
        fh = fopen(fn_in, "rb");

        if(fh == NULL) {
                fprintf(stderr, "FATAL: %s : cannot open binmat, returning NULL\n", __func__ );
                return NULL;
        }

        fread(&magicFront, sizeof(uint64_t), 1, fh);


        if(magicFront != MAGIC_HEAD) {
                fprintf(stderr, "FATAL: %s : binmat is corrupt, returning NULL\n", __func__ );
                return NULL;
        }

        struct matrix * m = malloc(sizeof(struct matrix));

        fread(&m->n1, sizeof(datum), 1, fh);
        fread(&m->n2, sizeof(datum), 1, fh);

        m->dat = malloc(sizeof(fwf*)*m->n1);

        datum index = 0;

        for(; index < m->n1; index += 1)
        {
                m->dat[index] = malloc(sizeof(fwf)*m->n2);
                fread(m->dat[index], sizeof(fwf), m->n2, fh);
        }
        fread(&magicBack, sizeof(uint64_t), 1, fh);

        if(magicBack != MAGIC_TAIL) {
                fprintf(stderr, "FATAL: %s : binmat is corrupt, returning NULL\n", __func__ );
                destroy_matrix(m); m = NULL; return m;
        }
        return m;
}
