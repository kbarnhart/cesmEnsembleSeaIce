#!/bin/bash
div=6
for i in `seq 1 85`;
do
	j=$i
	k=$((j/div))
	cmd2="qsub -l nodes=compute-1-$k:ppn=1 runDLM_$i.sh"
    echo $cmd2
    $cmd2
done 




