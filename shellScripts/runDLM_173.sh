#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1801,2000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
