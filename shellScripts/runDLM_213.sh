#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9801,10000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 