#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9201,9400, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
