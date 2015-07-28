#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(501,1000, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
