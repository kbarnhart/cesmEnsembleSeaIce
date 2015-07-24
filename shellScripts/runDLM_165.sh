#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(201,400, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
