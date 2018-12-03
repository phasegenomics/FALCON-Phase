# FALCON-Phase  ![Build Status](https://travis-ci.com/phasegenomics/FALCON-Phase.svg?token=qJUQGbDRUX3LsN3id6ky&branch=master) [![Analytics](https://ga-beacon.appspot.com/UA-90096627-2/FALCON-Phase/blob/master/README.md)](https://github.com/igrigorik/ga-beacon)

FALCON-Phase integrates [PacBio](https://www.pacb.com/) long-read assemblies with [Phase Genomics](https://phasegenomics.com/) Hi-C data to create phased, diploid, chromosome-scale scaffolds.

![FP logo](https://github.com/phasegenomics/FALCON-Phase/blob/master/logo/FP.png)

Read about the method and its performance in our [preprint](http://biorxiv.org/cgi/content/short/327064v1).

Version 2 of this method will be maintained in the pb-assembly bioconda package. See [here](https://github.com/PacificBiosciences/pb-assembly) for details of the full suite of tools.

## Recommended Usage Case :monkey:

FALCON-Phase was developed to solve the problem of haplotype switching in diploid genome assemblies. It has been tested on mammalian assemblies and performs well on an F1 bull with [0.7% heterozygosity](https://www.biorxiv.org/content/early/2018/02/26/271486). It has not been thoroughly tested on samples with very low or very high heterozygosity. In principle, it can be applied to outbred diploid organisms with less than 5% divergence between maternal and paternal haplotypes; when run on organisms with higher heterozygsity determining homology between haplotypes is ineffective as currently implemented. For samples with low heterozygosity, you should determine how much of your genome has "unzipped" by comparing the lengths of your haplotigs to those of your primary contigs. If the total haplotig length is small, FALCON-Phase will have limited utility.

To run the pipeline you need a [FALCON-Unzip](https://github.com/PacificBiosciences/pb-assembly#general-overview) assembly and Hi-C data. See [PacBio](https://www.pacb.com/calculator-whole-genome-sequencing/) recommendations for assembly coverage. For Hi-C we suggest 100 million read pairs per 1 Gb of genome length, with adjustments for genome complexity and library quality.

FALCON-Phase can be used to phase haplotype blocks within a contig and contigs within scaffolds. See our [preprint](http://biorxiv.org/cgi/content/short/327064v1) and SFAF [slides](http://pb-falcon.readthedocs.io/en/latest/_downloads/FALCON-Phase_SBKingan_23May2018SFAF.pdf) for details.


## Dependencies :rage2:

We did our best to minimize dependencies, but there are a number of standard bioinformatics tools required by the pipeline.
The version numbers of the dependencies are listed below, but new/older versions should work, but are untested. The required binaries are specified in the config.json file.  
+ **Python      (3.6)**                 -  Running Snakemake
+ **NumPy       (1.14.2)**              -  Pipeline python scripts
+ **Snakemake   (3.6)**                 -  Running the pipeline interactively or on a cluster 
+ **BWA         (v0.7.17)**             -  Mapping the Hi-C to the minced contigs
+ **Mummer 4    (4.0.0)**               -  Mapping the haplotigs (h contigs)
+ **BEDTools    (2.27.1)**              -  Creating AB pair index 
+ **HSTLIB      (1.5 or greater)**      -  Internal dependency (bundled with FALCON-Phase)
+ **SAMTOOLS    (1.5 or greater)**      -  Indexing fasta files



## Installation :floppy_disk:

Cloning the repository downloads the pipeline, the source code, and HTSLIB. 

`git clone --recursive https://github.com/phasegenomics/FALCON-Phase.git ; cd FALCON-Phase/src ; make`

After running this command you should see the executable binary `FALCON-Phase/bin/falcon-phase`.

## Running the test dataset :horse_racing:

We have provided a small pseudo test dataset ([sample info](https://www.ncbi.nlm.nih.gov/assembly/GCA_002021895.1/)) to get the pipeline running. 


1. Install the pipeline (as shown in the Install Process)

2. Clean the headers of the test FALCON-Unzip assembly files and generate the name_mapping.txt file. We have provides a perl script to do this for you: `scrub_names.pl` in the bin directory. The script removes the "|arrow" suffix on the FASTA headers and makes new files with the ".cleaned" suffix on the base name.

Usage:

```scrub_names.pl p-contigs.fa h-contigs.fa > name_mapping.txt```

The name_mapping.txt file maps each haplotig to a primary contig and is required to for the mummer steps and to build the snakemake DAG. It looks like this:

```000000F	000000F_001
000000F	000000F_002
000000F	000000F_003
000001F	000001F_001
000001F	000001F_002
000001F	000001F_003
```

3. Edit the config files

If you are running on a cluster that uses modules, you can load the binaries with a series of commands in the `config.sh` file. Below is an example `config.sh`. 

```
module load snakemake
module load bwa/0.7.17
module load bedtools/2.27.1
module load samtools/1.7
module load mummer/4.0.0
```

Next, edit the `config.json` in the pipeline folder, filling out the paths to the dependencies, and sample information. The paths to the required binaries in the config.json can be determined using `which` (e.g. `which samtools` returns the path to the binary).

The tables below explains the fields in the config file. 

### Enviromental setup :vhs:

| Key           | Value                                             | Explanation     |
| ------------- |:-------------------------------------------------:|:---------------|
| env           | config.sh                                         | This file is sourced by the shell and loads enviromental variables |
| CPU           | 2                                                 |  default number of CPUs to use |
| nucmer        | /path/to/mummer/bin/nucmer                        |  Path to nucmer |
| delta-filter  | /path/to/mummer/bin/delta-filter                  |  Path to delta-filter |
| show-coords   | /path/to/mummer/bin/show-coords                   |  Path to show-coords |
| samtools      | /path/to/samtools                                 |  Path to samtools |
| hp            | /path/to/FALCON-Phase/bin/coords2hp.py            |  Path to coords2hp.py |
| hpfilt        | /path/to/FALCON-Phase/bin/filt_hp.py              |  Path to filt_hp.py |
| falcon_phase  | /path/to/FALCON-Phase/bin/falcon-phase            |  Path to falcon-phase | 
| falcon_oi     | /path/to/FALCON-Phase/bin/primary_contig_index.pl |  Path to primary_contig_index.pl | 
| bedtools      | /path/to/bedtools                                 |  Path to bedtools | 
| bwa           | /path/to/bwa                                      |  Path to bwa |
| bwa           | cpu: 24                                           |  number of CPUs for bwa mapping |


### Sample setup :link:

| Key          | Value                                 | Explanation                                      |
| ------------ |:-------------------------------------:|:------------------------------------------------ |
| name         | test                                  | The name of the sample, most output files will have this prefix |
| output_format| pseudohap                             | The desired output format ("unzip" or "pseudohap") |
| min_aln_len  | 3000                                  | The minimal alignment length to consider during haplotig placement |
| p_to_h       | /path/to/name_mapping_tiny.txt        | A file that maps the haplotig names to primary contig names |
| p_ctgs       | /path/to/cns_p_ctg.clean.fasta | Path to CLEANED primary contigs |
| h_ctgs       | /path/to/cns_h_ctg.clean.fasta | Path to CLEANED haplotigs |
| r1           | /path/to/FALCON-Phase/test_dataset/S3HiC_R1.fastq  | Hi-C read-pair 1 |
| r2           | /path/to/FALCON-Phase/test_dataset/S3HiC_R2.fastq  | Hi-C read-pair 2 |
| enzyme       | GATC                                   | The restriction enzyme recognition seq used for Hi-C library prep |
| iter         | 10000                                  | The number of iterations for phasing algorithm, 10e7 is recommended |


4. Phew! Take a deep breath, you're almost done.

5a. Run the `snakemake` pipeline locally. To do so, on my computer, I need to be in a python 3.6 virtual environment with `NumPy` libraries loaded. This differs between machines, so you'll need to get that sorted out. Here's the command to run the FALCON-Phase snakemake pipeline:

```
(py36)$ snakemake -s snakefile --verbose -p
```

Snakemake will start printing information to your screen about the process. This example dataset only takes a few minutes to run, but a "real" job will take longer and should be run in the background. Once everthing is done you should see four files in the `output` direstory.

5b. Run the `snakemake` on a cluster. Snakemake can be run on a cluster and submit jobs the scheduler. We have included an SGE cluster config file as an example. The table below explains the fields in the `cluster.config.sge.json`. 

| Field   | Key     | Value           | Explanation                                |
| ------- | ------- |:---------------:|:------------------------------------------ |
| default | S       | /bin/bash       | interpreting shell |
|         | N       | falcon-phase    | job name |
|         | P       | -cwd            | directory |
|         | Q       | default         | queue name |
|         | CPU     | -pe smp 2       | default number of cores requested |
|         | E       | qsub_log/       | dir for stderr |
|         | O       | qsub_log/       | dir for stdout |
| aln     | CPU     | -pe smp 24      | number of cores requested for bwa mapping |

__NOTES__: the number of CPUs specified in the `cluster.config.sge.json` should match that in environmetal settings in config.json. Make sure you create the `qsub_log` dir before launching the job!!

Below is the command to run snakemake on PacBio's SGE cluster. This command runs 50 concurrent jobs and pulls the other `qsub` parameters from the `cluster.config.sge.json` file.

```
snakemake -j 50 --cluster-config cluster.config.sge.json --cluster "qsub -S {cluster.S} -N {cluster.N} {cluster.P} -q {cluster.Q} {cluster.CPU} -e {cluster.E} -o {cluster.O} -V" -s snakefile --verbose -p --latency-wait 60
```


## Directory Structure and Pipeline Output :cow2:

The FALCON-Phase pipeline has a number of steps that are implemented by rules in the `snakefile`. `Snakemake` implements the rules from bottom to top in the `snakefile`. Details about each step can be found in our [preprint](http://biorxiv.org/cgi/content/short/327064v1) and are summarized in the figure below:

![FP Workflow](https://github.com/phasegenomics/FALCON-Phase/blob/master/logo/Fig1_sbk.pdf)

### Starting Assembly (Workflow Step 1)

FALCON-Phase is only compatible with assemblies generated by FALCON-Unzip. Primary contigs have names like 000123F and haplotigs have names like 000123F_001.

### Haplotig Placement (Workflow Step 2)

The haplotig placement file specifies where each haplotig aligns to its primary contig and is used to define phase blocks in the FALCON-Unzip assembly. Making the haplotig placement file involved several steps. First, `nucmer` is used to align all FALCON-Unzip haplotigs to their primary contig. This step produces delta files for each primary contig.

        job_dir/
        ├── delta_files/
            ├── test.000000F.delta/             # delta file of haplotigs aligned to primary contig 000000F
            ├── test.000001F.delta/             # delta file of haplotigs aligned to primary contig 000001F

Delta files are filtered in the next step using `delta-filter`:

        job_dir/
        ├── filtered_delta_files/
            ├── test.000000F.delta.filt/        # filtered delta file of haplotigs aligned to primary contig 000000F
            ├── test.000001F.delta.filt/        # filtered delta file of haplotigs aligned to primary contig 000001F

Alignment coordinate files are made with `show-coords`:

        job_dir/
        ├── coords_files/
            ├── test.000000F.coords/        # coords file for primary contig 000000F
            ├── test.000001F.coords/        # coords file for primary contig 000001F

The haplotig placement file is made by calling two scripts, `coords2hp.py` then `filt_hp.py`. The file specifying the phase block pairing is contained in the haplotig_placement directory and is produced in the mincing stage (see below).

        job_dir/
        ├── haplotig_placement_file/
            ├── test.hbird.hp.txt/             # unfiltered haplotig placement file
            ├── test.hbird.filt_hp.txt/        # final haplotig placement file
            ├── test.AB_pairs.txt/             # pairing of A-B haplotigs (phase blocks)

### Mincing (Workflow Step 3)

Once the haplotig placement file and A-B phase block pairings are done, the primary contigs are minced at phase block boundaries. Mincing allows Hi-C reads to be mapped to each pair of phase block so that the density of Hi-C connections can be used to assign blocks to the same phase.

        job_dir/
        ├── mince/
            ├── test.A_haplotigs.bed/    # BED file for A minced haplotigs (original FALCON-Unzip haplotigs)
            ├── test.B_haplotigs.bed/    # BED file for B minced haplotigs (corresponding phase blocks on FALCON-Unzip primary contigs)
            ├── test.collapsed_haplotypes.bed/  #  BED file for collapsed haplotigs (non-Unzipped regions of primary contigs)
            ├── test.A_haplotigs.fasta/  # FASTA file for A minced haplotigs 
            ├── test.B_haplotigs.fasta/  # FASTA file for B minced haplotigs
            ├── test.collapsed_haplotypes.fasta/  #  FASTA file for collapsed haplotypes
            ├── test.minced.fasta/        # Concatenated FASTA (in this order: A_haplotigs, B_haplotigs, collapsed)
            ├── test.BC.bed/              # sorted BED of B and collapsed
            ├── B_haplotigs_merged.bed/   # merging of overlapping B haplotigs (used to define collapsed regions)


### Mapping of Hi-C Reads (Workflow Step 4)

Hi-C reads are mapped to the minced FASTA file using `bwa mem`, streamed and processed in `samtools` and filtered with a FALCON-Phase utility. This is usually the most time consuming step and typically runs overnight on a mammalian genome with 24-48 cores.

        job_dir/
        ├── hic_mapping/
                ├── test.unfiltered.bam/     # Hi-C reads mapped to minced FASTA
                ├── test.filtered.bam/       # filtered Hi-C reads mapped to minced FASTA
        

### Phasing (Workflow Step 5)

The binary matrix and index file are input to the FALCON-Phase algorithm to assign phase to each A-B phase block pair. The index file is created by the `primary_contig_index.pl` script.

        job_dir/
        ├── phasing/
                ├── test.binmat/            # binary matric of normalized Hi-C mappings
                ├── test.ov_index.txt/      # index file of ordering of minced contigs on each primary contig and A-B phase block pairings
                ├── test.results.txt/        # phase assignment for each A-B phase block pair
	        ├── test.seq.txt/        # summary of RE cut sites for each minced fragment


The `results.txt` file has the following columns:

1. primary contig ID
2. phase 0 haplotig ID
3. phase 1 haplotig ID
4. proportion of simulations in this configuration
5. normalized count of reads mapped to phase 0 haplotig
6. normalized count of reads mapped to phase 1 haplotig
7. numeric index in fasta for sequence 0
8. numeric index in fasta for sequence 1


### Emission of Phased Haplotigs (Workflow Step 6)

You have two options for the final output of FALCON-Phase, which is generated by the `emit_haplotigs.pl` script. `pseudohap` output is two phased full-length pseudo-haplotypes for each primary contig. `unzip` style output is a primary contig file and haplotigs file with the phase-switch errors corrected. "Phase 0" output corresponds to the primary contigs and "Phase 1" to the haplotigs.

        job_dir/
        ├── output/
                ├── test_data.phased.0.bed/    # list of minced contigs used to reconstruct phase 0 output
                ├── test_data.phased.1.bed/    # list of minced contigs used to reconstruct phase 1 output
                ├── test_data.phased.0.fasta/   # phase 0 output
                ├── test_data.phased.1.fasta/   # phase 1 output

For the `pseudohap` the final output FASTA headers look like:

```
>000000F_0
>000001F_0
>000000F_1
>000001F_1
```

Phase is specified by the `_0` or `_1` suffix on the original FALCON-Unzip primary contig IDs.

For the `unzip` style output, the fasta headers are the same as what you are used to with FALCON-Unzip:

Primary contigs:
```
>000000F
>000001F
```
Haplotigs:
```
>000000F_001
>000000F_003
>000000F_002
>000001F_003
>000001F_001
>000001F_002
```





