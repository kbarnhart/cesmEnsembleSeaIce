#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1401,1600, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 