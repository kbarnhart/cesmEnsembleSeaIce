#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9501,10000, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
