#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2801,3000, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
