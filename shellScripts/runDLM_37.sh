#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(7201,7400, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
