#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7501,8000, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
