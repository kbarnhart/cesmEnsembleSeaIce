#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(4401,4600, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
