#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(6401,6600, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
