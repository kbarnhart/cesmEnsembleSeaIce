#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9501,10000, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
