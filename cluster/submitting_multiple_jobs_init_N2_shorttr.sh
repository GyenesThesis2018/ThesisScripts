#!/bin/bash

for As in `seq 1 6`; do
	qsub -v a=$As -N sig_proc_job_gs_$As submit_init_and_subset_shorttr_N2.pbs
done