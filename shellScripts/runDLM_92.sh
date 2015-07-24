#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7401,7600, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
