#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7201,7400, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
