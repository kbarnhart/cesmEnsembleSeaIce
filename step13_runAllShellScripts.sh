#!/bin/bash

for i in `seq 1 88`;
do
	cmd2="qsub -l nodes=1 -l walltime=48:00:00 runDLM_$i.sh"
    echo $cmd2
    $cmd2
done