## run locally
#snakemake -s snakefile --verbose -p

## run distributed on sge cluster
## 20 jobs submitted at a time
## 60s wait time
## -V option ports env to compute nodes
#snakemake -j 20 --cluster-config cluster.config.sge.json --cluster "qsub -S {cluster.S} -N {cluster.N} {cluster.P} -q {cluster.Q} {cluster.CPU} -e {cluster.E} -o {cluster.O} -V" -s snakefile --verbose -p --latency-wait 60
