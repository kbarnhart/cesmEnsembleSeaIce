#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7001,7500, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
