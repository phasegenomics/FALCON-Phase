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


#include "count_motif.h"

int print_sequenceInfo(FILE * stream, struct sequenceInfo * sinfo)
{
	if(sinfo == NULL) return 1;

	int i = 0;
	int j = 0;

	fprintf(stream, "%s\t%s\t%s\n", "#sequence", "length", "cutsites"  );
	for(; i < sinfo->nseq; i++) {
		fprintf(stream, "%s\t%i", sinfo->dat[i].name, sinfo->dat[i].slen);
		for(j = 0; j < sinfo->ncutters; j++) {
			fprintf(stream, "\t%s:%i", sinfo->cutter_motifs[j], sinfo->dat[i].cutcount[j] );
		}
		fprintf(stream, "\n");
	}
	return 0;
}

int32_t count_cutsite( char * motif, faidx_t * fai, const char * seqName, int seq_idx)
{


	int mlen = strlen(motif);


	int i = 0;
	while(i < strlen(motif)) {
		motif[i] = toupper(motif[i]);
		i++;
	}

	int n = faidx_nseq(fai);

	int32_t motif_n = 0;

	char * sname = seqName;

	if(seqName == NULL) {
		if(seq_idx >= n) {
			fprintf(stderr, "FATAL: sequence index > number of sequences in fasta : in : %s \n ", __func__ );
			free (sname);
			return -1;
		}
		sname = faidx_iseq(fai, seq_idx);
	}

	if(faidx_has_seq(fai, sname) == 0) {
		free (sname);
		return -1;
	}

	int slen = faidx_seq_len(fai, sname);

	int jnk = 0;

	char * seq = faidx_fetch_seq(fai, sname, 0, slen-1, &jnk);

	i = 0;
	int len =  strlen(seq);

	while(i < len) {
		seq[i] = toupper(seq[i]);
		i++;
	}
	char * tmp_ptr_for  = seq;
	char * tmp_ptr_last = seq;

	while(1) {
		if(tmp_ptr_for == NULL ) break;
		tmp_ptr_for = strstr(tmp_ptr_last, motif);

		tmp_ptr_last = mlen + tmp_ptr_for;


		if(tmp_ptr_for == NULL ) break;
		motif_n++;

	}

	free(seq);

	return motif_n;
}


int count_cutsite_runner(char ** argv, int argc)
{

	fprintf(stderr, "INFO: running: %s on : %s\n", __func__, argv[2] );


	faidx_t * fai = fai_load(argv[2]);
	if(fai == NULL) {
		fprintf(stderr, "FATAL: fasta: %s could not be loaded/indexed in %s\n", argv[2], __func__  );
		exit(1);
	}

	int nseq = faidx_nseq(fai);
	int32_t motifCount = 0;
	int i    = 0;
	int nmof = 3;

	fprintf(stdout, "##seqid");
	for(; nmof < argc; nmof++) {
		fprintf(stdout, "\t%s", argv[nmof]);
	}
	fprintf(stdout, "\n");

	for(; i < nseq; i++) {
		nmof = 3;
		fprintf(stdout, "%s", faidx_iseq(fai, i));
		for(; nmof < argc; nmof++) {

			motifCount   = 0;
			motifCount  += count_cutsite(argv[nmof], fai, NULL, i);

			if(motifCount < 0) {
				fprintf(stderr, "FATAL something went wrong in %s %i\n", __func__, i );
				fai_destroy(fai);
				return 1;
			}
			fprintf(stdout, "\t%i", motifCount );
		}
		fprintf(stdout, "\n");
	}


	fai_destroy(fai);
	return 0;
}



struct sequenceInfo * load_seq_info(char * fasta, char * motif){

	// si is the number of enzymes determined by comma count
	int i, si;
	si = 1;

	for(i = 0; i < strlen(motif); i++) {
		if(motif[i] == ',') si++;
		if(motif[i] == 'N') return NULL;
		if(motif[i] == 'n') return NULL;

	}

	fprintf(stderr, "INFO: loading sequence information\n" );

	struct sequenceInfo * seqInfo;
	seqInfo = malloc(sizeof(struct sequenceInfo));

	seqInfo->ncutters      = si;
	seqInfo->cutter_motifs = malloc(sizeof(char * )*si);

	int s = 0;
	int e = 0;

	for(i = 0; i < si; i++) {
		seqInfo->cutter_motifs[i] = get_next_word(motif, &s, &e, ',');
		fprintf(stderr, "INFO: added motif: %s\n", seqInfo->cutter_motifs[i] );
	}

	seqInfo->fai = fai_load(fasta);
	if(seqInfo->fai == NULL) {
		fprintf(stderr, "FATAL: fasta: %s could not be loaded/indexed in %s\n", fasta, __func__  );
		free(seqInfo);
		exit(1);
	}

	int nseq = faidx_nseq(seqInfo->fai);

	seqInfo->dat = malloc(sizeof(struct sequence)*nseq);
	seqInfo->nseq = nseq;

	for(i = 0; i < nseq; i++) {

		seqInfo->dat[i].cutcount = malloc(sizeof(uint32_t)*seqInfo->ncutters);

		seqInfo->dat[i].name     = faidx_iseq(seqInfo->fai, i);
		seqInfo->dat[i].slen     = faidx_seq_len(seqInfo->fai, seqInfo->dat[i].name);

		int j = 0;
		for(j = 0; j < seqInfo->ncutters; j++) {
			seqInfo->dat[i].cutcount[j] = count_cutsite(seqInfo->cutter_motifs[j], seqInfo->fai, NULL, i);
		}
	}

	fprintf(stderr, "INFO: loaded sequence information\n" );
	return seqInfo;
}

struct sequenceInfo * destroy_sequence_info(struct sequenceInfo * seqInfo){

	if(seqInfo == NULL) return NULL;

	int i = 0;
	for(i = 0; i < seqInfo->ncutters; i++) {
		free(seqInfo->cutter_motifs[i]);
	}

	for(i = 0; i < seqInfo->nseq; i++) {
		free(seqInfo->dat[i].cutcount);
	}

	fai_destroy(seqInfo->fai);


	seqInfo = NULL;
	return seqInfo;
}
