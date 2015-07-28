#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2001,2500, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
