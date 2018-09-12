.PHONY: all clean

OBJ = matrix.o bam_filt.o count_motif.o dna.o run_phasing.o string_parser.o 
SRC = matrix.c bam_filt.c count_motif.c dna.c run_phasing.c string_parser.c

CFLAG=-Wall -O3

INCLUDE=-I../htslib -Isrc 

all: ../bin/falcon-phase
	@printf "\n\n[INFO] build OK, found falcon-phase!\n\n"

../htslib/libhts.a:
	cd ../htslib && $(MAKE)

../bin/falcon-phase: ../bin ../htslib/libhts.a main.c $(OBJ)
	gcc $(CFLAG) main.c $(OBJ) $(INCLUDE) ../htslib/libhts.a  -lbz2 -llzma -lz -lm -lpthread  -o ../bin/falcon-phase

$(OBJ) : %.o : %.c ../htslib/libhts.a
	gcc $(CFLAG) $(INCLUDE) -c $< -o $@


clean:
	@rm -rf *.o; rm -rf ../bin/falcon-phase
