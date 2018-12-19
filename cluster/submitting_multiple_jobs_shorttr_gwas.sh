#!/bin/bash

for As in `seq 106 106`; do for Bs in `seq 2 2`; do
	qsub -v a=$As,b=$Bs -N sig_proc_job_$As submit_array_shorttr_gwas.pbs
done
done