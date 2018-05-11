#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "matrix.h"
#include "bam_filt.h"
#include "run_phasing.h"

int flag_counts[6] = {0,0,0,0,0,0};

void print_usage(void){
								fprintf(stderr, "\nusage: falcon-phase <command> [options] \n\n");
								fprintf(stderr, "\ncommands:\n - bam2 - converts alignments to several useful hi-c formats.\n");
								fprintf(stderr, "   + usage: falcon-phase bam2 [binmat|lachesis|juicer|counts] input output\n");
								fprintf(stderr, "   + details: \n");
								fprintf(stderr, "      The input file format is automatically determined [cram|bam|sam].\n");
								fprintf(stderr, "      The output is written to the fineame provided, no extention.\n");
								fprintf(stderr, "\n        \n - bamfilt - filter a hi-c bam.\n");
								fprintf(stderr, "   + usage: falcon-phase bamfilt input.[cram|bam|sam] output.bam\n");
								fprintf(stderr, "\n        \n - cutsites - count cutsites per seqid.\n");
								fprintf(stderr, "   + usage: falcon-phase cutsites input.fasta ATGC TGCA ...\n");

}

int process_pair_juicer( bam_hdr_t *header, FILE * fh){

								char * rname1 = bam_get_qname(rp.read1);
								char * rname2 = bam_get_qname(rp.read2);

								if(strcmp(rname1, rname2) != 0) return 0;

								if(rp.read1->core.qual == 0 ) return 2;
								if(rp.read2->core.qual == 0 ) return 2;

								uint8_t *sa = bam_aux_get(rp.read1, "SA");
								uint8_t *xa = bam_aux_get(rp.read1, "XA");
								if(xa != 0 || sa != 0) return 2;

								sa = bam_aux_get(rp.read2, "SA");
								xa = bam_aux_get(rp.read2, "XA");
								if(xa != 0 || sa != 0) return 2;

								int strandA = 0;
								int strandB = 0;

								if((rp.read1->core.flag & 16)) {
																strandA = 16;
								}
								if((rp.read2->core.flag & 16)) {
																strandB = 16;
								}

								char* tName = header->target_name[rp.read1->core.tid];
								char* qName = header->target_name[rp.read2->core.tid];
								if(rp.read1->core.tid <= rp.read2->core.tid) {
																fprintf(fh, "%i %s %i %i %i %s %i %i 1 - - 1  - - -\n", strandA, tName, \
																								rp.read1->core.pos, 0, strandB, qName, rp.read1->core.mpos, 1);
								}
								else{
																fprintf(fh, "%i %s %i %i %i %s %i %i 1 - - 1  - - -\n", strandB, qName, \
																								rp.read1->core.mpos, 1, strandA, tName, rp.read1->core.pos, 0);
								}

								return 1;
}

int process_pair_lachesis( bam_hdr_t *header, struct matrix * lp){


								char * rname1 = bam_get_qname(rp.read1);
								char * rname2 = bam_get_qname(rp.read2);

								// mates must have the same read neame
								if(strcmp(rname1, rname2) != 0)
								{
																fprintf(stderr, "FATAL read pairs out of sync: %s != %s\n", rname1, rname2);
																return -1;
								}



								if(add_link(lp, rp.read1->core.tid, rp.read2->core.tid) < 1) {
																fprintf(stderr, "FATAL: add_link\n");
																exit(1);
								}
								if(add_link(lp, rp.read2->core.tid, rp.read1->core.tid) < 1) {
																fprintf(stderr, "FATAL: add_link\n");
																exit(1);
								}

								return PASS;
}

int bam2lachesis(const char * fn_in, const char * fn_out, int type){

								fprintf(stderr, "INFO: converting bam to lachesis on %s\n", fn_in );

								rp.read1 = bam_init1();
								rp.read2 = bam_init1();

								htsFile * h = hts_open(fn_in, "r");

								const htsFormat *fmt = hts_get_format(h);

								int hts_close(htsFile *fp);

								fprintf(stderr, "INFO: detected %s filetype\n", hts_format_file_extension(fmt));

								samFile *in = 0;
								bam_hdr_t *header = NULL;

								if((in = sam_open_format(fn_in, "r", fmt)) == 0) {
																printf("FATAL: failed to open \"%s\" for reading\n", fn_in);
																return 1;
								}

								if ((header = sam_hdr_read(in)) == 0) {
																fprintf(stderr, "FATAL: failed to read the header \"%s\" \n", fn_in);
																return 1;
								}

								uint32_t nSeqs = header->n_targets;

								struct matrix * lpc = init_matrix(nSeqs, nSeqs);

								fprintf(stderr, "INFO: reading file \"%s\"\n", fn_in);

								int r1;
								int r2;
								long int c = 0;

								while (1) {

																c+=1;

																r1 = sam_read1(in, header, rp.read1);
																r2 = sam_read1(in, header, rp.read2);

																if( r1 < 0 || r2 < 0 ) break;

																int flag = 0;

																flag = process_pair_lachesis(header, lpc);


																if(flag < 0) {
																								fprintf(stderr, "FATAL: something went wrong in process_pair\n");
																								exit(1);
																}
																if((c % 1000000) == 0) {
																								fprintf(stderr, "INFO: parsed %ld read pairs\n", c);
																}
								}

								if(type == 1) {
																print_matrix(lpc, fn_out);
								}
								if(type == 0) {
																freeze_matrix(lpc, fn_out);
								}
								if(type == 2) {
																FILE * fh;
																fh = fopen(fn_out, "wb");
																if(fh == NULL) return 1;

																fprintf(stderr, "INFO: matrix_size: %i by %i\n", lpc->n1, lpc->n2);

																datum ol = 0;
																datum il = 0;

																for(; ol < lpc->n1; ol++) {
																								il = 0;
																								for(; il < lpc->n2; il++) {
																																if(lpc->dat[ol][il] != 0) {
																																								fprintf(fh, "%s\t%s\t%f\n",
																																																header->target_name[ol],
																																																header->target_name[il],
																																																lpc->dat[ol][il]);
																																}
																								}
																}
								}

								destroy_matrix(lpc);

								hts_close(in);
								bam_destroy1(rp.read1);
								bam_destroy1(rp.read2);
								bam_hdr_destroy(header);


								return 0;

}

int bam2juicer(const char * fn_in, const char * fn_out){

								FILE * fh;
								fh = fopen(fn_out, "wb");
								if(fh == NULL) return 0;

								fprintf(stderr, "INFO: converting bam to juicer on %s\n", fn_in );

								rp.read1 = bam_init1();
								rp.read2 = bam_init1();

								htsFile * h = hts_open(fn_in, "r");

								const htsFormat *fmt = hts_get_format(h);

								int hts_close(htsFile *fp);

								fprintf(stderr, "INFO: detected %s filetype\n", hts_format_file_extension(fmt));

								samFile *in = 0;
								bam_hdr_t *header = NULL;

								if((in = sam_open_format(fn_in, "r", fmt)) == 0) {
																printf("FATAL: failed to open \"%s\" for reading\n", fn_in);
								}
								if ((header = sam_hdr_read(in)) == 0) {
																fprintf(stderr, "FATAL: failed to read the header \"%s\" \n", fn_in);
								}


								fprintf(stderr, "INFO: reading file \"%s\"\n", fn_in);

								int r1;
								long int c = 0;

								while (1) {

																c+=1;

																r1 = sam_read1(in, header, rp.read1);
																r1 = sam_read1(in, header, rp.read2);
																//        fprintf(stderr, "INFO: r1 %i read pairs\n", r1);

																if( r1 < 0 ) break;
																if(process_pair_juicer(header, fh) < 1) {
																								fprintf(stderr, "FATAL: something went wrong in process_pair\n");
																								exit(1);
																}
																if((c % 1000000) == 0) {
																								fprintf(stderr, "INFO: parsed %ld read pairs\n", c);
																}
								}

								hts_close(in);
								bam_destroy1(rp.read1);
								bam_destroy1(rp.read2);
								bam_hdr_destroy(header);

								return 0;
}

int main(int argc, char **argv){

								if(argc < 2) {
																print_usage();
																return 1;
								}

								if(strcmp(argv[1], "bin2d") == 0) return run_binner(argv, argc);
								if(strcmp(argv[1], "cutsites") == 0) return count_cutsite_runner(argv, argc);
								if(strcmp(argv[1], "readcount") == 0) return count_reads(argv[2]);
								if(strcmp(argv[1], "bam2") == 0) {
																if(strcmp(argv[2], "binmat")  == 0) return bam2lachesis(argv[3], argv[4], 0);
																if(strcmp(argv[2], "counts")  == 0) return bam2lachesis(argv[3], argv[4], 2);
																if(strcmp(argv[2], "lachesis")== 0) return bam2lachesis(argv[3], argv[4], 1);
																if(strcmp(argv[2], "juicer")  == 0) return bam2juicer(argv[3],      argv[4]);
								}
								if(strcmp(argv[1], "phase") == 0) return run_phasing(argv, argc);
								if(strcmp(argv[1], "bamfilt") == 0) return filter_bam(argv, argc);


								fprintf(stderr, "FATAL unknown command %s\n", argv[2]);
								print_usage();
								return 1;
}
