#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9401,9600, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
