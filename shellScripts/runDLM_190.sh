#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(5201,5400, 'SIdata.sh.RCP85.mat'), exit" -nodesktop -nosplash 
