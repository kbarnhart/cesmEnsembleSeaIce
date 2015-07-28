#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3501,4000, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
