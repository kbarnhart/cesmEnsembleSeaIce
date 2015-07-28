#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(5001,5500, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
