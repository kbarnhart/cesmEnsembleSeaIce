#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2501,3000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
