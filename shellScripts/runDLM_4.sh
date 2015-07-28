#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(1501,2000, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
