#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9601,9800, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
