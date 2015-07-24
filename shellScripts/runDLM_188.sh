#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4801,5000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
