#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(8001,8200, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
