#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7601,7800, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
