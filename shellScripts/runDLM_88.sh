#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10501,10914, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
