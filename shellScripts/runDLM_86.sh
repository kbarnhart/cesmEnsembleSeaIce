#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(6201,6400, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
