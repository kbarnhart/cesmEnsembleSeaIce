#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10001,10200, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
