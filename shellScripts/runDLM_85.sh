#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9001,9500, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
