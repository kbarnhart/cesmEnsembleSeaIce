#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10201,10400, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
