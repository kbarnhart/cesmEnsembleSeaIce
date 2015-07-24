#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2601,2800, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
