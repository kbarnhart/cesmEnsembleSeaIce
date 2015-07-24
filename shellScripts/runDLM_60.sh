#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1001,1200, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
