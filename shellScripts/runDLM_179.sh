#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3001,3200, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
