#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2001,2200, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 
