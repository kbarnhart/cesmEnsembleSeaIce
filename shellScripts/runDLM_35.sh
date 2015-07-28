#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(6001,6500, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
