#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10001,10500, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
