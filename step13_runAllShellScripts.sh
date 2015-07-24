#!/bin/bash
div=6
for i in `seq 1 218`;
do
	j=$i
	k=$((j/div))
	cmd2="qsub -l nodes=1 runDLM_$i.sh"
    echo $cmd2
    $cmd2
done 




