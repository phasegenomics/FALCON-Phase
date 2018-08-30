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


#include "run_phasing.h"
#include "string_parser.h"

void print_phasing_usage(void)
{
	fprintf(stderr, "\n\nusage: falcon-phase phase [options] -f your.fasta -b your.binmat -m GATC -p sample \n\n");
	fprintf(stderr, "   + required:\n");
	fprintf(stderr, "      -f fasta file.\n");
	fprintf(stderr, "      -b binmat file.\n");
	fprintf(stderr, "      -m comma separated list of motifs, eg. GATC or GATC,TTTC. No degenerate bases! (list all possible motifs)\n");
	fprintf(stderr, "      -p output prefix\n");
	fprintf(stderr, "      -i overlap index file\n");
	fprintf(stderr, "   + options: \n");
	fprintf(stderr, "      -n number of iterations [2,000,000]\n");
	fprintf(stderr, "      -u number of  burnin iterations [1,000,000]\n");
	fprintf(stderr, "      -s seed value for rand() [time]\n");
	fprintf(stderr, "      -v dump large data [false]\n\n");
}

/**
 * parse_command_line_phasing description : load globally scoped options struct
 * @param  argv - argv from main
 * @param  argc - argc from main
 * @return      [description]
 */
int parse_command_line_phasing(char ** argv, int argc)
{

	global_opts_phasing.fasta      = NULL;
	global_opts_phasing.binmat     = NULL;
	global_opts_phasing.nsweeps    = 2000000;
	global_opts_phasing.burnin     = 1000000;
	global_opts_phasing.motif      = NULL;
	global_opts_phasing.prefix     = NULL;
	global_opts_phasing.index      = NULL;
	global_opts_phasing.verbose    = 0;
	global_opts_phasing.seed       = -1;


	int c;
	const char    * short_opt = "hf:b:n:m:p:v:i:u:s:";
	struct option long_opt[] =
	{
		{"help",          no_argument,       NULL, 'h'},
		{"index",         required_argument, NULL, 'i'},
		{"fasta",         required_argument, NULL, 'f'},
		{"bin",           required_argument, NULL, 'b'},
		{"niters",         optional_argument, NULL,'n'},
		{"cutoff",        optional_argument, NULL, 'u'},
		{"motif",         required_argument, NULL, 'm'},
		{"prefix",        required_argument, NULL, 'p'},
		{"verbose",       no_argument,       NULL, 'v'},
		{"seed",          optional_argument, NULL, 's'},
		{NULL,            0,                 NULL, 0  }
	};

	while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
	{
		switch(c)
		{
		case 'v':
		{
			global_opts_phasing.verbose = 1;
			break;
		}
		case 'h':
		{
			print_phasing_usage();
			exit(1);
		}
		case 'm':
		{
			global_opts_phasing.motif = optarg;
			break;
		}
		case 'p':
		{
			global_opts_phasing.prefix = optarg;
			break;
		}
		case 'f':
		{
			global_opts_phasing.fasta = optarg;
			break;
		}
		case 'i':
		{
			global_opts_phasing.index = optarg;
			break;
		}
		case 'd':
		{
			global_opts_phasing.damping = atof(optarg);
			break;
		}
		case 'b':
		{
			global_opts_phasing.binmat = optarg;
			break;
		}
		case 'n':
		{
			global_opts_phasing.nsweeps = atoi(optarg);
			break;
		}
		case 'u':
		{
			global_opts_phasing.burnin = atoi(optarg);
			break;
		}
		case 's':
		{
			global_opts_phasing.seed = atoi(optarg);
			break;
		}


		default:
		{
			/* invalid option */
			fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",
			        argv[0], optopt);
			break;
		}
		}
	}

	if(global_opts_phasing.burnin > global_opts_phasing.nsweeps) {
		fprintf(stderr, "FATAL: burn in must be less than the number of sweeps \n");
		print_phasing_usage();
		exit(1);
	}

	fprintf(stderr, "INFO: RUNNING: ");
	int i = 0;

	for(; i<argc; i++)
	{
		fprintf(stderr, "%s ", argv[i]);
	}
	fprintf(stderr, "\n");
	return 0;
}

/*
 * normalize_matrix_phasing description: takes the matrix.h object and normalizes the hi-c counts based on the number of restriction sites on each sequence.
 * @param  mat - matrix object pointer defined in matrix.h
 * @param  seq - sequenceInfo object pointer defined in count_motif.c. It contains the sequence name/length/cutsite counts
 * @return     [description]
 */
int normalize_matrix_phasing(struct matrix * mat, struct sequenceInfo * seq)
{

	double * mins = (double *)malloc(sizeof(double)*mat->n1);

	datum i,j,k;

	for(i = 0; i < mat->n1; i++) {
		mins[i] = DBL_MAX;
	}

	int n_cutters = seq->ncutters;

	for(i = 0; i < mat->n1; i++) {
		for(j = 0; j < mat->n2; j++) {
			if(mat->dat[i][j] == 0) {
				continue;
			}
			double cutsite_sum = 0;
			for(k = 0; k < n_cutters; k++) {
				cutsite_sum += seq->dat[i].cutcount[k];
				cutsite_sum += seq->dat[j].cutcount[k];
			}

			if(cutsite_sum == 0) {
				mat->dat[i][j] = 0;
				fprintf(stderr, "WARNING: seq: %s and seq: %s have zero cutsites, setting hi-c counts to 0\n", seq->dat[i].name, seq->dat[j].name );
				continue;
			}

			mat->dat[i][j] = mat->dat[i][j] / (double)cutsite_sum;

		}
	}

	return 0;
}

/**
 * parse_index_file description: this ugly function parses the overlap index file using a custom string parser defined in
 * @param  fn [description]
 * @return    [description]
 */
struct line_data * parse_index_file(const char * fn){

	struct line_data * head = NULL;
	struct line_data * last = NULL;

	FILE * fp;
	fp = fopen(fn, "r");

	fprintf(stderr, "INFO: parsing index file %s\n", fn);

	if(fp == NULL) {
		fprintf(stderr, "FATAL: did not load index file");
		return 0;
	}

	size_t len = 0;
	ssize_t read;
	char * line = NULL;

	while( (read = getline(&line, &len, fp)) != -1) {

		int s = 0;
		int e = 0;

		char * group        = get_next_word(line, &s, &e, 9);
		char * master_index = get_next_word(line, &s, &e, 9);
		char * overlap      = get_next_word(line, &s, &e, 9);
		overlap[strlen(overlap)-1] = '\0';

		int nctg = 1;
		int i;
		for(i = 0; i < strlen(master_index); i++) {
			if(master_index[i] == ',') nctg++;
		}

		int noverlap = 1;
		for(i =0; i < strlen(overlap); i++) {
			if(overlap[i] == ',') noverlap++;
		}



		struct line_data * ld = (struct line_data *) malloc(sizeof(struct line_data));
		ld->group    = group;
		ld->indices  = (int *) malloc(sizeof(int) * nctg);
		ld->overlaps = (struct ctg_overlap * ) malloc(sizeof(struct ctg_overlap) * noverlap);
		ld->has_overlap = 0;
		if(noverlap >= 1) ld->has_overlap = 1;
		ld->n_overlap = noverlap;
		ld->n_ctg     = nctg;
		ld->before    = NULL;
		ld->after     = NULL;

		// parse indicies
		s = 0;
		e = 0;
		i = 0;
		char * indstr = NULL;
		while((indstr = get_next_word(master_index, &s, &e, ',')) != NULL) {
			ld->indices[i] = atoi(indstr);
			free(indstr);
			i++;
		}

		if(strcmp("NA", overlap) == 0) ld->has_overlap = 0;

		if(ld->has_overlap) {
			s = 0;
			e = 0;
			i = 0;
			indstr = NULL;
			while((indstr = get_next_word(overlap, &s, &e, ',')) != NULL) {
				int ss = 0;
				int ee = 0;

				char * indbstr = get_next_word(indstr, &ss, &ee, ':');


				ld->overlaps[i].first = atoi(indbstr);
				free(indbstr);

				indbstr = get_next_word(indstr, &ss, &ee, ':');
				ld->overlaps[i].second = atoi(indbstr);
				free(indbstr);


				int first = rand() & 1;

				ld->overlaps[i].label_count_one    = 0;
				ld->overlaps[i].label_count_two    = 0;
				ld->overlaps[i].label_one = first;
				ld->overlaps[i].label_two = first ^ 1;


				assert(ld->overlaps[i].label_one !=   ld->overlaps[i].label_two);
				i++;
			}
		}

		if(head == NULL) {
			head = ld;
		}
		else{
			ld->before  = last;
			last->after =   ld;
		}
		last = ld;

		free(master_index);
		free(overlap);
	}

	free(line);
	fclose(fp);

	return head;

}

/*
 * phase_likelihood description: Calculates the within vs between hi-c signal for all past diced haplotigs, given the past configuration
 * @param  ov - overlap index
 * @param  mv - normalized hi-c matrix
 * @param  i  - index of current overlap (index)
 * @return    [description]
 */
double phase_likelihood(struct ctg_overlap * ov,
                        struct matrix * mv,
                        int i)
{

	// we set the within and between values to avoid a division by zero
	double within  = 0.000001;
	double between = 0.000001;

	int j;
	for(j = i -1; j >= 0; j--) {

		if(ov[i].label_one == ov[j].label_one) {
			within  += mv->dat[ov[i].first][ov[j].first];
			within  += mv->dat[ov[i].second][ov[j].second];

			between += mv->dat[ov[i].first][ov[j].second];
			between += mv->dat[ov[i].second][ov[j].first];
		}
		else{
			between  += mv->dat[ov[i].first][ov[j].first];
			between  += mv->dat[ov[i].second][ov[j].second];

			within += mv->dat[ov[i].first][ov[j].second];
			within += mv->dat[ov[i].second][ov[j].first];
		}
	}

	return within / (within + between);
}

/**
 * print_local: a function to print data for troubleshooting. I don't remeber what this even does so don't mess with it.
 * @param  ov [description]
 * @param  mv [description]
 * @param  n  [description]
 * @return    [description]
 */
int _print_local(struct ctg_overlap * ov,
                 struct matrix *mv,
                 int n){
	int i,j;

	fprintf(stdout, "graph overlaps {\nrankdir=\"BT\";\n");

	for(i = 0; i < n; i++) {
		fprintf(stdout, "%i -> %i [label=\"%.3f\",weight=\"%.3f\",style=dotted,color=blue];\n", ov[i].first,
		        ov[i].second, mv->dat[ov[i].first][ov[i].second], mv->dat[ov[i].first][ov[i].first] +mv->dat[ov[i].second][ov[i].second]);

		for(j = i + 1; j < n; j++) {
			fprintf(stdout, "%i -- %i [label=\"%.3f\"weight=\"%.3f\"];\n", ov[i].first,
			        ov[j].first, mv->dat[ov[i].first][ov[j].first], mv->dat[ov[i].first][ov[j].first]);

			fprintf(stdout, "%i -- %i [label=\"%.3f\",weight=\"%.3f\"];\n", ov[i].first,
			        ov[j].second, mv->dat[ov[i].first][ov[j].second], mv->dat[ov[i].first][ov[j].second]);


			fprintf(stdout, "%i -- %i [label=\"%.3f\",weight=\"%.3f\"];\n", ov[i].second,
			        ov[j].first, mv->dat[ov[i].second][ov[j].first], mv->dat[ov[i].second][ov[j].first]);


			fprintf(stdout, "%i -- %i [label=\"%.3f\",weight=\"%.3f\"];\n", ov[i].second,
			        ov[j].second, mv->dat[ov[i].second][ov[j].second], mv->dat[ov[i].second][ov[j].second]);

		}
	}
	fprintf(stdout, "{rank = same; ");
	for(i = 0; i < n; i++) {
		fprintf(stdout, " %i;", ov[i].first);
	}
	fprintf(stdout, "}\n");
	fprintf(stdout, "{rank = same; ");
	for(i = 0; i < n; i++) {
		fprintf(stdout, " %i;", ov[i].second);
	}
	fprintf(stdout, "}\n");
	fprintf(stdout, "}\n");
	return 0;
}

/**
 * stochastic_phasing description : this is the function that does the phasing. It loops over the overlap index structure per primary contig and stochastically searches for the best haplotig configuration based on the hi-c data.
 * @param  ld  - pointer to line_data array, each element is a haplotype.
 * @param  si  - pointer to the sequenceInfo structure. see count_motif.h
 * @param  lc  - pointer to hi-c matrix. see matrix.h
 * @param  n   - number of primary contigs that have haplotigs.
 * @param  res - FILE * to print the results to.
 * @return     - zero upon success
 */
int stochastic_phasing(struct line_data * ld,
                       struct sequenceInfo * si,
                       struct matrix * lc, int n,
                       FILE * results)
{

	struct line_data * tmpld = ld;

	while(tmpld != NULL) {
		fprintf(stderr, "INFO: working on group %s\n", tmpld->group );
		if(!tmpld->has_overlap) {
			tmpld = tmpld->after;
			continue;
		}

		uint32_t iter, j;

		for(iter = 0; iter < global_opts_phasing.nsweeps; iter++) {

			for(j = 0; j < tmpld->n_overlap; j++) {
				double r = (double)rand() / (double)RAND_MAX;
				double prob =  phase_likelihood(tmpld->overlaps, lc, j);

				if(r > prob && j != 0) {
					tmpld->overlaps[j].label_one ^= 1;
					tmpld->overlaps[j].label_two ^= 1;
				}
				assert(tmpld->overlaps[j].label_one != tmpld->overlaps[j].label_two);

				// if we've gone past the burning we start counting the orientation of each overlap (phase block)
				if(iter > global_opts_phasing.burnin  && tmpld->overlaps[j].label_one == 1) tmpld->overlaps[j].label_count_one += 1;
				if(iter > global_opts_phasing.burnin  && tmpld->overlaps[j].label_one == 0) tmpld->overlaps[j].label_count_two += 1;
			}
		}

		for(j = 0; j < tmpld->n_overlap; j++) {
			const char * s1 = faidx_iseq(si->fai, tmpld->overlaps[j].first);
			const char * s2 = faidx_iseq(si->fai, tmpld->overlaps[j].second);

			double v1 = lc->dat[tmpld->overlaps[j].first][tmpld->overlaps[j].first];
			double v2 = lc->dat[tmpld->overlaps[j].second][tmpld->overlaps[j].second];

			double p = tmpld->overlaps[j].label_count_one / (double)(tmpld->overlaps[j].label_count_one + tmpld->overlaps[j].label_count_two);
			if(p <= 0.5) {
				fprintf(results, "%s %s %s %f %.4f %.4f %i %i\n", tmpld->group, s2, s1, 1-p, v2, v1, tmpld->overlaps[j].second, tmpld->overlaps[j].first );
			}
			else{
				fprintf(results, "%s %s %s %f %.4f %.4f %i %i\n", tmpld->group, s1, s2, p, v1, v2, tmpld->overlaps[j].first, tmpld->overlaps[j].second);
			}
		}
		assert(tmpld != tmpld->after);
		tmpld = tmpld->after;
	}
	return 0;
}

/*
 * run_phasing description : main
 * @param  argv [description]
 * @param  argc [description]
 * @return      [description]
 */
int run_phasing(char ** argv, int argc)
{



	int parse_flag = parse_command_line_phasing(argv, argc);
	if(parse_flag > 0) return 1;

	if(global_opts_phasing.fasta  == NULL ||
	   global_opts_phasing.binmat == NULL ||
	   global_opts_phasing.motif  == NULL ||
	   global_opts_phasing.prefix == NULL   ) {
		fprintf(stderr, "FATAL: missing one of the input files %s %s %s %s %s\n", global_opts_phasing.fasta,  global_opts_phasing.binmat, global_opts_phasing.motif,  global_opts_phasing.prefix,  global_opts_phasing.index);
		print_phasing_usage();
		return 1;
	}

	int seed = time(NULL);
	if(global_opts_phasing.seed  != -1) {
		seed = global_opts_phasing.seed;
	}
	srand(seed);


	char * cutsite_fn = malloc(strlen(global_opts_phasing.prefix) + 12);
	char * data_fn    = malloc(strlen(global_opts_phasing.prefix) + 12);
	char * results_fn = malloc(strlen(global_opts_phasing.prefix) + 12);

	strcat(cutsite_fn, global_opts_phasing.prefix); strcat(cutsite_fn, ".seqs.txt");
	strcat(data_fn,    global_opts_phasing.prefix); strcat(data_fn,    ".data.txt");
	strcat(results_fn, global_opts_phasing.prefix); strcat(results_fn, ".results.txt");

	FILE * cutsites;
	FILE * results;

	cutsites = fopen(cutsite_fn, "w");
	results  = fopen(results_fn, "w");

	// read the index file
	struct line_data * index_data = parse_index_file(global_opts_phasing.index);

	struct line_data * tp = index_data;

	int n = 0;
	while(tp != NULL) {
		tp = tp->after;
		n += 1;
	}

	// load the sequence info
	struct sequenceInfo * sInfo = load_seq_info(global_opts_phasing.fasta, global_opts_phasing.motif);

	if(sInfo == NULL) {
		fprintf(stderr, "FATAL: problem counting cutsites in FASTA\n");
		return(1);
	}

	if(sInfo == NULL) return 1;
	print_sequenceInfo(cutsites, sInfo);

	// read the matrix
	struct matrix * link_counts = thaw_matrix(global_opts_phasing.binmat);

	// normlize the hi-c data
	normalize_matrix_phasing(link_counts, sInfo);

	// the wrapper function to run the phasing algorthim
	stochastic_phasing(index_data, sInfo, link_counts, n, results);

	fclose(cutsites);
	fclose(results);

	// cleanup allocated structures
	destroy_matrix(link_counts);
	destroy_sequence_info(sInfo);

	return 0;

}
