#!/bin/bash 
cd /home/barnhark/seaIceEmergence 
matlab -r "seaIceEmergeDLM_v1(2601,2800, 'SIdata.sh.RCP45.mat'), exit" -nodesktop -nosplash 