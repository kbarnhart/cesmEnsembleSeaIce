#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4001,4200, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
