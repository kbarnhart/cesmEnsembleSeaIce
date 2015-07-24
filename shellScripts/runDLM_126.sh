#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3401,3600, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
