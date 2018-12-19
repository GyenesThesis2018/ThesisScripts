#!/bin/bash

for As in `seq 79 79`; do for Bs in `seq 5 5`; do
	qsub -v a=$As,b=$Bs -N sig_proc_job_$As submit_array_shorttr_N2.pbs
done
done