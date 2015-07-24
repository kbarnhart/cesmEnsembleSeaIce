#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(5601,5800, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
