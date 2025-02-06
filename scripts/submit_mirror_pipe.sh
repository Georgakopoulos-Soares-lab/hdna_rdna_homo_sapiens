#!/bin/bash

snakemake --snakefile submit_mirror_repeat.smk \
         --rerun-triggers mtime \
         --keep-incomplete \
         --rerun-incomplete \
         --latency-wait 45 \
         --keep-going \
         --cluster-config cluster_settings.yaml \
         --cluster "sbatch -p {cluster.partition} -t {cluster.time} --mem={cluster.mem} -c {cluster.ncpus} --nodes={cluster.nodes} -o jobOut/{cluster.jobName}-%j.out -J {cluster.jobName} -e jobOut/{cluster.jobName}-%j.err" \
         -j 12
