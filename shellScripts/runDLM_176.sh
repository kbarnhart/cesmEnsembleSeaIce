#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2401,2600, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
