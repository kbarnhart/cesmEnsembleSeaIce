#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(3201,3400, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 