#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1001,1500, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
