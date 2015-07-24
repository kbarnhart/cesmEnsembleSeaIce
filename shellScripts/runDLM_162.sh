#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10601,10800, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
