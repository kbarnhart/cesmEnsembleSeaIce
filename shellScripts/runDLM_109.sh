#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1,200, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
