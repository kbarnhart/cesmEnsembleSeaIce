#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(6601,6800, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
