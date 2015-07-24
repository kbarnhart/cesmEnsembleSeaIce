#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1201,1400, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
