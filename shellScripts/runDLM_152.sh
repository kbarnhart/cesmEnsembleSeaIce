#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(8601,8800, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
