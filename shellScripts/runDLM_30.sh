#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(5801,6000, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
