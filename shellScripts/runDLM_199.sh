#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7001,7200, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
