#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(10601,10726, 'SIdata.nh.RCP85.mat'), exit" -nodesktop -nosplash 
