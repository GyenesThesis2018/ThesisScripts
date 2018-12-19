#!/bin/bash

for As in `seq 2 2`; do
	qsub -v a=$As -N sig_proc_job_gs_$As submit_init_and_subset_shorttr_gwas.pbs
done