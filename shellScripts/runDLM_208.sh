#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(8801,9000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 