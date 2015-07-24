#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(5401,5600, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
