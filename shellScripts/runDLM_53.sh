#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4001,4500, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
