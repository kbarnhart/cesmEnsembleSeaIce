#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(8201,8400, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
