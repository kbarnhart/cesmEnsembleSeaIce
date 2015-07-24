#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3201,3400, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
