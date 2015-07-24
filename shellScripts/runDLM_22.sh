#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4201,4400, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
