#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(6801,7000, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
