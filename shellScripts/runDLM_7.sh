#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3001,3500, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
