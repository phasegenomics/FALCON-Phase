#include "bam_filt.h"


struct bam_filt_opts {
	char * in;
	char * out;
	char * exclude_str;
	int * exclude_lookup;
	int * tid_remap;
	int * covered;
	int include;
	int mq_filter;
	int ed_filter;
	int bin_filter;
	uint32_t min_t_len;
} bamfilt_global_opts;


void print_bam_filt_usage(void){
	fprintf(stderr, "\n    usage: matlock bamfilt [options] -i input.[cram|bam|sam] -o output.bam \n\n"              );
	fprintf(stderr, "         Required:\n"                                                                             );
	fprintf(stderr, "    -i -         Input file.\n"                                                        );
	fprintf(stderr, "    -o -         Ouput file.\n\n"                                                      );
	fprintf(stderr, "         Options:\n"                                                     );
	fprintf(stderr, "    -h -          Print help statement.\n"                                                   );
	fprintf(stderr, "    -m - <INT>    MapQ filter. [20]    \n"                                                             );
	fprintf(stderr, "    -e - <INT>    Max edit distance. [5]\n"                                                          );
	fprintf(stderr, "    -l - <INT>    Min target seq-length. [0]\n"                                                          );
	fprintf(stderr, "    -x - <STRING> Comma separated list of seqids to exclude/include. [exclude]\n"                                   );
	fprintf(stderr, "                  This option should be used with the binary flag 64 (-f 64).\n"                                             );
	fprintf(stderr, "    -y -          incude -x rather than exclude [exclude]\n"                                                    );
	fprintf(stderr, "    -f - <INT>    Binary flag filter:\n\n"                                                    );
	fprintf(stderr, "                  SAME_SEQID  =  2\n"                                                                 );
	fprintf(stderr, "                  LOW_MAPQ    =  4\n"                                                                 );
	fprintf(stderr, "                  XA_SA       =  8\n"                                                                 );
	fprintf(stderr, "                  NM          = 16\n"                                                                 );
	fprintf(stderr, "                  SMALLCONTIG = 32\n"                                                                 );
	fprintf(stderr, "                  EXCLUDE     = 64\n"                                                                 );
	fprintf(stderr, "                  SA_ONLY     = 128\n"                                                                );
	fprintf(stderr, "                  UNMAPPED    = 256\n"                                                                );
	fprintf(stderr, "                  DUPLICATE   = 1024\n"                                                                );
	fprintf(stderr, "                  The default is 1300 =  LOW_MAPQ | NM | DUPLICATE | UNMAPPED \n\n"                                            );

}

char *inputString(FILE* fp, size_t size){
//The size is extended by the input with the value of the provisional
	char *str;
	int ch;
	size_t len = 0;
	str = realloc(NULL, sizeof(char)*size);//size is start size
	if(!str) return str;
	while(EOF!=(ch=fgetc(fp)) && ch != '\n') {
		str[len++]=ch;
		if(len==size) {
			str = realloc(str, sizeof(char)*(size+=16000));
			if(!str) return str;
		}
	}
	str[len++]='\0';

	return realloc(str, sizeof(char)*len);
}


int parse_length_exclude(bam_hdr_t * header){
	int i = 0;
	for(; i < header->n_targets; i++) {
		if(header->target_len[i] < bamfilt_global_opts.min_t_len) {
			bamfilt_global_opts.exclude_lookup[i] = 1;
		}
	}
	return 0;
}

int parse_target_list(bam_hdr_t *header){

	bamfilt_global_opts.exclude_lookup = malloc(sizeof(int) * header->n_targets);
	memset(bamfilt_global_opts.exclude_lookup, 0, header->n_targets * sizeof(int));
	if(bamfilt_global_opts.include == 1) {

		int i = 0;
		for(; i < header->n_targets; i++) {
			bamfilt_global_opts.exclude_lookup[i] = 1;
		}
	}

	if(header->n_targets == 0) {
		fprintf(stderr, "WARNING : no sequences in bam header : in : %s\n", __func__);
	}

	// strtok will cause a seg fault if you giv it a NULL at least for Unbuntu
	if(bamfilt_global_opts.exclude_str == NULL) return 0;

	char * pch = strtok (bamfilt_global_opts.exclude_str,",");
	uint32_t nseq = 0;
	while (pch != NULL)
	{
		nseq = 0;
		for(; nseq < header->n_targets; nseq++) {
			if(strcmp(pch, header->target_name[nseq]) == 0) {
				if(bamfilt_global_opts.include == 1) {
					//      fprintf(stderr, "INFO: including %s : in : %s \n", header->target_name[nseq], __func__ );
					bamfilt_global_opts.exclude_lookup[nseq] = 0;
				}
				else{
					//      fprintf(stderr, "INFO: excluding %s : in : %s \n", header->target_name[nseq], __func__ );
					bamfilt_global_opts.exclude_lookup[nseq] = 1;
				}
			}
		}
		pch = strtok (NULL, ",");
	}
	return 0;
}

int parse_bamfilt_command_line(char ** argv, int argc)
{

	bamfilt_global_opts.in             = NULL;
	bamfilt_global_opts.out            = NULL;
	bamfilt_global_opts.mq_filter      = 20;
	bamfilt_global_opts.ed_filter      = 5;
	bamfilt_global_opts.bin_filter    |= LOW_MQ | NM | DUPLICATE | UNMAPPED;
	bamfilt_global_opts.min_t_len      = 0;
	bamfilt_global_opts.exclude_str    = NULL;
	bamfilt_global_opts.exclude_lookup = NULL;
	bamfilt_global_opts.include        = 0;



	int c;
	const char    * short_opt = "hm:e:f:l:x:i:o:y";
	struct option long_opt[] =
	{
		{"help",          no_argument,       NULL, 'h'},
		{"include",       no_argument,       NULL, 'y'},
		{"mapq",          optional_argument, NULL, 'm'},
		{"input",         required_argument, NULL, 'o'},
		{"output",        required_argument, NULL, 'i'},
		{"edit",          optional_argument, NULL, 'e'},
		{"flag",          optional_argument, NULL, 'f'},
		{"length",        optional_argument, NULL, 'l'},
		{"exclude",       optional_argument, NULL, 'x'},
		{NULL,            0,                 NULL, 0  }
	};

	while((c = getopt_long(argc, argv, short_opt, long_opt, NULL)) != -1)
	{
		switch(c)
		{
		case 'h':
		{
			print_bam_filt_usage();
			exit(1);
		}
		case 'm':
		{
			bamfilt_global_opts.mq_filter = atoi(optarg);
			bamfilt_global_opts.bin_filter |= LOW_MQ;
			break;
		}
		case 'y':
		{
			bamfilt_global_opts.include = 1;
			break;
		}
		case 'e':
		{
			bamfilt_global_opts.ed_filter = atoi(optarg);
			bamfilt_global_opts.bin_filter |= NM;
			break;
		}
		case 'f':
		{
			bamfilt_global_opts.bin_filter = atoi(optarg);
			break;
		}
		case 'l':
		{
			bamfilt_global_opts.min_t_len = atoi(optarg);
			bamfilt_global_opts.bin_filter |= SMALLCONTIG;
			break;
		}
		case 'x':
		{
			bamfilt_global_opts.exclude_str = optarg;
			if( access( bamfilt_global_opts.exclude_str, F_OK ) != -1 ) {
				FILE * exclude_fn = fopen(bamfilt_global_opts.exclude_str, "r");
				bamfilt_global_opts.exclude_str = bamfilt_global_opts.exclude_str = inputString(exclude_fn, 1000);
				close(exclude_fn);
				bamfilt_global_opts.bin_filter |= EXCLUDE;
			}
			break;
		}
		case 'i':
		{

			bamfilt_global_opts.in = optarg;
			break;
		}
		case 'o':
		{
			bamfilt_global_opts.out = optarg;

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

	fprintf(stderr, "INFO: RUNNING: ");
	int i = 0;

	for(; i<argc; i++)
	{
		fprintf(stderr, "%s ", argv[i]);
	}
	fprintf(stderr, "\n");
	return 0;
}

int _process_pair_bamfilt(bam_hdr_t *header, int mapqf, int editf, uint32_t tlen_min){

	int binflag = 0;

	if(rp.read1->core.flag & 4 || rp.read2->core.flag & 4  ) {
		return UNMAPPED;
	}

	if(bamfilt_global_opts.exclude_lookup[rp.read1->core.tid] == 1) binflag |= EXCLUDE;
	if(bamfilt_global_opts.exclude_lookup[rp.read2->core.tid] == 1) binflag |= EXCLUDE;

	if(header->target_len[rp.read1->core.tid] < tlen_min) binflag |= SMALLCONTIG;
	if(header->target_len[rp.read2->core.tid] < tlen_min) binflag |= SMALLCONTIG;

	// mates must map to different contigs
	if(rp.read1->core.tid == rp.read2->core.tid) binflag |= SAME_SEQID;

	// alignment quality must be greater than mapqf;
	if(rp.read1->core.qual < mapqf ) binflag |= LOW_MQ;
	if(rp.read2->core.qual < mapqf ) binflag |= LOW_MQ;

	char * rname1 = bam_get_qname(rp.read1);
	char * rname2 = bam_get_qname(rp.read2);

	// mates must have the same read neame
	if(strcmp(rname1, rname2) != 0)
	{
		fprintf(stderr, "FATAL read pairs out of sync: %s != %s\n", rname1, rname2);
		exit(1);
	}

	// multimapping bad boys
	uint8_t *sa = bam_aux_get(rp.read1, "SA");
	uint8_t *xa = bam_aux_get(rp.read1, "XA");
	if(xa != 0 || sa != 0) binflag |= XA_SA;
	if(sa != 0) binflag |= SA_ONLY;

	sa = bam_aux_get(rp.read2, "SA");
	xa = bam_aux_get(rp.read2, "XA");
	if(xa != 0 || sa != 0) binflag |= XA_SA;
	if(sa != 0) binflag |= SA_ONLY;


	uint8_t *nm1 = bam_aux_get(rp.read1, "NM");
	uint8_t *nm2 = bam_aux_get(rp.read2, "NM");

	if(nm1 != NULL && nm2 != NULL) {

		if(bam_aux2i(nm1) > editf) binflag |= NM;
		if(bam_aux2i(nm2) > editf) binflag |= NM;
	}

	// filter out duplicates
	if (rp.read1->core.flag & DUPLICATE) binflag |= DUPLICATE;
	if (rp.read2->core.flag & DUPLICATE) binflag |= DUPLICATE;

	return binflag;
}


int filter_bam(char ** argv, int argc){

	parse_bamfilt_command_line(argv, argc);


	if(bamfilt_global_opts.in == NULL || bamfilt_global_opts.out == NULL) {

		fprintf(stderr, "FATAL: missing -i:%s or -o:%s\n", bamfilt_global_opts.in, bamfilt_global_opts.out );
		print_bam_filt_usage();
		exit(1);
	}

	int flag_counts[5000];
	memset(&flag_counts, 0, 5000* sizeof(int));

	htsFile * h = hts_open(bamfilt_global_opts.in, "r");

	const htsFormat *fmt = hts_get_format(h);

	int hts_close(htsFile *fp);

	fprintf(stderr, "INFO: detected %s filetype\n", hts_format_file_extension(fmt));

	samFile *in = 0;
	bam_hdr_t *header = NULL;

	if((in = sam_open_format(bamfilt_global_opts.in, "r", fmt)) == 0) {
		printf("FATAL: failed to open \"%s\" for reading\n", bamfilt_global_opts.in);
		return 1;
	}

	if ((header = sam_hdr_read(in)) == 0) {
		fprintf(stderr, "FATAL: failed to read the header \"%s\" \n", bamfilt_global_opts.in);
		return 1;
	}


	samFile *output = sam_open(bamfilt_global_opts.out, "wb");

	parse_target_list(header);
	/*
	   This must go after parse_target_list.
	 */
	parse_length_exclude(header);


	if(sam_hdr_write(output, header) != 0) {
		fprintf(stderr, "FATAL: failed to read the write \"%s\" \n", bamfilt_global_opts.out);
		return 1;
	}


	rp.read1 = bam_init1();
	rp.read2 = bam_init1();

	int okay;
	int r1;
	int r2;
	long int c = -1;

	char * rn1;
	char * rn2;

	while (1) {

		c+=1;

		r1 = sam_read1(in, header, rp.read1);
		r2 = sam_read1(in, header, rp.read2);


		if( r1 < 0 || r2 < 0 ) break;

		int flag = 0;

		flag = _process_pair_bamfilt(header,
		                             bamfilt_global_opts.mq_filter,
		                             bamfilt_global_opts.ed_filter,
		                             bamfilt_global_opts.min_t_len);

		flag_counts[flag & SAME_SEQID ]  += 1;
		flag_counts[flag & LOW_MQ     ]  += 1;
		flag_counts[flag & XA_SA      ]  += 1;
		flag_counts[flag & NM         ]  += 1;
		flag_counts[flag & SMALLCONTIG]  += 1;
		flag_counts[flag & EXCLUDE    ]  += 1;
		flag_counts[flag & SA_ONLY    ]  += 1;
		flag_counts[flag & DUPLICATE  ]  += 1;
		flag_counts[flag & UNMAPPED   ]  += 1;

		if((c % 1000000) == 0) {
			fprintf(stderr, "INFO: parsed %ld read pairs\n", c);
		}

		if((flag & bamfilt_global_opts.bin_filter) == 0) {
			flag_counts[200] += 1;

			okay = sam_write1(output, header, rp.read1);
			if(okay < -1) {
				fprintf(stderr, "[FATAL] issue writing bam in %s near line %s \n", __FILE__, __LINE__);
				exit(1);
			}
			okay = sam_write1(output, header, rp.read2);
			if(okay < -1) {
				fprintf(stderr, "[FATAL] issue writing bam in %s near line %s \n", __FILE__, __LINE__);
				exit(1);
			}
		}
	}

	bam_destroy1(rp.read1);
	bam_destroy1(rp.read2);
	bam_hdr_destroy(header);
	sam_close(output);
	free(bamfilt_global_opts.exclude_lookup);
	free(bamfilt_global_opts.tid_remap);
	free(bamfilt_global_opts.covered);


	fprintf(stderr, "STATS: mate pair that passed filtering:...... %i %f%%\n", flag_counts[200], (flag_counts[200] / (double)c) * 100);
	fprintf(stderr, "STATS: mate pair that are not mapped:........ %i %f%%\n", flag_counts[256], (flag_counts[256] / (double)c) * 100);

	fprintf(stderr, "STATS: mate pair with low mapq (%i):......... %i %f%%\n", bamfilt_global_opts.mq_filter, flag_counts[LOW_MQ], (flag_counts[LOW_MQ] / (double)c)*100 );
	fprintf(stderr, "STATS: mate pair with XA or SA tag:.......... %i %f%%\n", flag_counts[XA_SA], (flag_counts[XA_SA] / (double)c)*100 );
	fprintf(stderr, "STATS: mate pair with same seqid:............ %i %f%%\n", flag_counts[SAME_SEQID],(flag_counts[SAME_SEQID]/(double)c)*100);
	fprintf(stderr, "STATS: mate pair with NM > %i:................ %i %f%%\n", bamfilt_global_opts.ed_filter, flag_counts[NM], (flag_counts[NM]/(double)c)*100 );
	fprintf(stderr, "STATS: mate pair on target sequences < %i Bp:. %i %f%%\n", bamfilt_global_opts.min_t_len, flag_counts[SMALLCONTIG], (flag_counts[SMALLCONTIG]/(double)c)*100 );
	fprintf(stderr, "STATS: mate pair in exclude list:............ %i %f%%\n", flag_counts[EXCLUDE], (flag_counts[EXCLUDE]/(double)c)*100 );
	fprintf(stderr, "STATS: mate pair with only SA tag:........... %i %f%%\n", flag_counts[SA_ONLY], (flag_counts[SA_ONLY] / (double)c)*100 );
	fprintf(stderr, "STATS: mate pair with a DUPLICATE flag:...... %i %f%%\n", flag_counts[DUPLICATE], (flag_counts[DUPLICATE] / (double)c)*100 );

	fprintf(stderr, "STATS: Total mate pair:...................... %ld\n", c);
	fprintf(stderr, "STATS: Total reads:.......................... %ld\n", c*2);


	return 0;
}
