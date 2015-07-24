#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(8601,8800, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
