#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1601,1800, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
