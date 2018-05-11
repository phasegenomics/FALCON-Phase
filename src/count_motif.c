#include "count_motif.h"

int print_sequenceInfo(FILE * stream, struct sequenceInfo * sinfo)
{
  if(sinfo == NULL) return 1;

  int i = 0;

    fprintf(stream, "%s\t%s\t%s\n", "#sequence", "length", "cutsites"  );
  for(; i < sinfo->nseq; i++){
    fprintf(stream, "%s\t%i\t%i\n", sinfo->dat[i].name, sinfo->dat[i].slen, sinfo->dat[i].cutcount  );
  }
  return 0;
}

int32_t count_cutsite( char * motif, faidx_t * fai, const char * seqName, int seq_idx)
{


  int mlen = strlen(motif);


  int i = 0;
  while(i < strlen(motif)){
    motif[i] = toupper(motif[i]);
    i++;
  }

  int n = faidx_nseq(fai);

  int32_t motif_n = 0;

  char * sname = seqName;

  if(seqName == NULL){
      if(seq_idx >= n){
        fprintf(stderr, "FATAL: sequence index > number of sequences in fasta : in : %s \n ", __func__ );
        free (sname);
        return -1;
      }
      sname = faidx_iseq(fai, seq_idx);
  }

  if(faidx_has_seq(fai, sname) == 0){
    free (sname);
    return -1;
  }

  int slen = faidx_seq_len(fai, sname);

  int jnk = 0;

  char * seq = faidx_fetch_seq(fai, sname, 0, slen-1, &jnk) ;

  i = 0;
  int len =  strlen(seq);

  while(i < len){
    seq[i] = toupper(seq[i]);
    i++;
  }
  char * tmp_ptr_for  = seq;
  char * tmp_ptr_last = seq;

  while(1){
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
    if(fai == NULL){
      fprintf(stderr, "FATAL: fasta: %s could not be loaded/indexed in %s\n", argv[2], __func__  );
      exit(1);
    }

    int nseq = faidx_nseq(fai);
    int32_t motifCount = 0;
    int i    = 0;
    int nmof = 3;

    fprintf(stdout, "##seqid");
    for(; nmof < argc; nmof++){
        fprintf(stdout, "\t%s", argv[nmof]);
    }
    fprintf(stdout, "\n");

    for(; i < nseq; i++){
      nmof = 3;
      fprintf(stdout, "%s", faidx_iseq(fai, i));
      for(; nmof < argc; nmof++){

        motifCount   = 0;
        motifCount  += count_cutsite(argv[nmof], fai, NULL, i);

        if(motifCount < 0){
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

  fprintf(stderr, "INFO: loading sequence information\n" );

  struct sequenceInfo * seqInfo;
  seqInfo = (struct sequenceInfo *)malloc(sizeof(struct sequenceInfo));
  seqInfo->fai = fai_load(fasta);
  if(seqInfo->fai == NULL){
    fprintf(stderr, "FATAL: fasta: %s could not be loaded/indexed in %s\n", fasta, __func__  );
    free(seqInfo);
    exit(1);
  }

  int nseq = faidx_nseq(seqInfo->fai);

  seqInfo->dat = (struct sequence *)malloc(sizeof(struct sequence)*nseq);
  seqInfo->nseq = nseq;

  int i = 0;
  for(; i < nseq; i++){

    seqInfo->dat[i].name     = faidx_iseq(seqInfo->fai, i);
    seqInfo->dat[i].slen     = faidx_seq_len(seqInfo->fai, seqInfo->dat[i].name);
    seqInfo->dat[i].cutcount = count_cutsite(motif, seqInfo->fai, NULL, i);
  }

  fprintf(stderr, "INFO: loaded sequence information\n" );
  return seqInfo;
}

struct sequenceInfo * destroy_sequence_info(struct sequenceInfo * seqInfo){

  if(seqInfo == NULL) return NULL;

 fai_destroy(seqInfo->fai);


  seqInfo = NULL;
  return seqInfo;
}
