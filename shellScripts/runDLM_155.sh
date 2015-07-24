#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(9201,9400, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
