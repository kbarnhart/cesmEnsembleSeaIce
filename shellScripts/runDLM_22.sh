#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10501,10726, 'SIdata.nh.RCP45.mat'), exit" -nodesktop -nosplash 
