#!/bin/bash

for As in `seq 1 6`; do
	qsub -v a=$As -N combj_$As submit_combining.pbs
done