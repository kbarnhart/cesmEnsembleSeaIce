#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4501,5000, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
