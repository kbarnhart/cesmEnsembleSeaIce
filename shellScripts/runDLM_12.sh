#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2201,2400, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 