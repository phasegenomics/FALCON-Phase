# FALCON-Phase  ![Build Status](https://travis-ci.com/phasegenomics/FALCON-Phase.svg?token=qJUQGbDRUX3LsN3id6ky&branch=master) [![Analytics](https://ga-beacon.appspot.com/UA-90096627-2/FALCON-Phase/blob/master/README.md)](https://github.com/igrigorik/ga-beacon)

FALCON-Phase integrates PacBio long-read assemblies with Phase Genomics Hi-C data to create phased, diploid, chromosome-scale scaffolds.

![FP logo](https://github.com/phasegenomics/FALCON-Phase/blob/master/logo/FP.png)

## Dependencies :rage2:

We did our best to minimize dependances, but there are a number of standard bioinformatics tools required by the pipeline.
The version numbers of the dependances are listed below, but new/older versions should work, but are untested. The required binaries are specified in the config.json file.  
+ **Python    (3.6)**             -  Running Snakemake
+ **Snakemake (3.6)**             -  Running the pipeline interactively or on a cluster 
+ **BWA       (v0.7.17)**         -  Mapping the Hi-C to the minced contigs
+ **Mummer 4  (4.0.0)**           -  Mapping the haplotigs (h contigs)
+ **BEDTools  (2.27.1)**          -  Creating AB pair index 
+ **HSTLIB    (1.5 or greater)**  -  Internal dependency (bundled with FALCON-phase)
+ **SAMTOOLS  (1.5 or greater)**  -  Indexing fasta files

In addition, the **NumPy** library is required.


## Install process :floppy_disk:

Cloning the repository downloads the pipeline, the source code, and HTSLIB. 

`git clone --recursive https://github.com/skingan/FALCON-phase.git ; cd FALCON-Phase-example/src ; make`

After running this command you should see the executable binary `FALCON-Phase/bin/falcon-phase`.

## Running the test dataset :running:

We have provided a small pseudo test dataset ([sample info](https://www.ncbi.nlm.nih.gov/assembly/GCA_002021895.1/)) to get the pipeline running. 


1. Install the pipeline (as shown in the Install Process)
2. Go into the pipeline folder.

```cd pipeline```

3. Edit the config.json, filling out the paths to the dependencies, and sample information. The tables below explains the fields in the config file. 

### Enviromental setup :vhs:

| Key           | Value                                             | Explanation     |
| ------------- |:-------------------------------------------------:|:---------------|
| env           | config.sh                                         | This file is sourced by the shell and loads enviromental variables |
| CPU           | number of CPUs to use                             |   2 |
| nucmer        | /path/to/mummer/bin/nucmer                        |   Path to nucmer |
| delta-filter  | /path/to/mummer/bin/delta-filter                  |   Path to delta-filter |
| show-coords   | /path/to/mummer/bin/show-coords                   |   Path to show-coords |
| samtools      | /path/to/samtools                                 |   Path to samtools |
| hp            | /path/to/FALCON-Phase/bin/coords2hp.py            |  Path to coords2hp.py |
| hpfilt        | /path/to/FALCON-Phase/bin/filt_hp.py              |  Path to filt_hp.py |
| falcon_phase  | /path/to/FALCON-Phase/bin/falcon-phase            |  Path to falcon-phase | 
| falcon_oi     | /path/to/FALCON-Phase/bin/primary_contig_index.pl |  Path to primary_contig_index.pl | 
| bedtools      | /path/to/bedtools                                 |  Path to bedtools| 
| bwa           | path: /path/to/bwa                                |  Path to bwa|
| bwa           | cpu: 24                                           |  number of CPUs|

If you already have binaries for the dependencies via cluster modules you will want to use the `config.sh` to load them. Below is an example `config.sh`. The paths to the binaries in the config.json will need to match where the cluster loads the modules (e.g. `which samtools` returns the path to the binary).

Example of config.sh
```
module load snakemake
module load bwa/0.7.17
module load bedtools/2.25.0
module load samtools/1.7
module load mummer/4.0.0
```



### Sample setup

| Key          | Value                                 | Explanation                                      |
| ------------ |:-------------------------------------:|:------------------------------------------------ |
| name         | test                                  | The name of the sample, most output files will have this prefix |
| min_aln_len  | 3000                                  | The minimal alignment length to consider during haplotig placement |
| h_to_p       | /path/to/name_mapping_tiny.txt        | A file that maps the haplotig names to primary contig names |
| p_ctgs       | /path/to/FALCON-Phase/test_dataset/cns_p_ctg.fasta | Path to primary contigs |
| h_ctgs       | /path/to/FALCON-Phase/test_dataset/cns_h_ctg.fasta | Path to haplotigs |
| r1           | /path/to/FALCON-Phase/test_dataset/S3HiC_R1.fastq  | Hi-C read-pair 1 |
| r2           | /path/to/FALCON-Phase/test_dataset/S3HiC_R2.fastq  | Hi-C read-pair 2 |
| enzyme       | GATC                                   | The restriction enzyme used for Hi-C library prep |
| iter         | 10000                                  | The number of iterations for phasing algorithm, 10e7 is recommended |


4. Phew! Take a deep breath, you're almost done.

5. Run the snakemake pipeline. To do so, on my computer, i need to be in a python 36 virtual environment. This differs between machines, so you'll need to get that sorted out. Here's the command to run the FALCON-Phase snakemake pipeline:

```
(py36) Zevs-MBP-2:pipeline zev$ snakemake -p
```

Snakemake will start printing information to your screen about the process. This example dataset only takes a few minutes to run, but a "real" job will take longer and should be run in the background. Once everthing is done you should see the final file `hbird.diploid_phased.fasta` in the pipeline folder.

Snakemake can be run on a cluster and submit jobs the scheduler. We have included an SGE cluster config file as an example. However, you'll need to setup the cluster config to match your cluster's settings. Below is the command to run snakemake on PacBio's SGE cluster.

```
snakemake -j 50 --cluster-config pb_cluster.config.sge.json --cluster "qsub -S {cluster.S} -N {cluster.N} {cluster.P} -q {cluster.Q} {cluster.CPU} -e {cluster.E} -o {cluster.O} -V" -s snakefile --verbose -p --latency-wait 60
```
