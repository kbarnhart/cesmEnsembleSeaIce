#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10801,10914, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
