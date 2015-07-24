#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4601,4800, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
