% fin='/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/justNSIF_ensembleAndBG.nc';
% nSIF_ensemble = ncread(fin,'nSIF_ensemble') ; 
% nSIF_background1850 = ncread(fin,'nSIF_background1850') ; 
% time_ensemble = ncread(fin,'time_ensemble') ; 
% time_background = ncread(fin,'time_background') ; 

load SIdata

NI=size(nSIF_ensemble,1);
NJ=size(nSIF_ensemble,2);
NM=size(nSIF_ensemble,4);
NY=size(nSIF_ensemble,3);
NBY=size(nSIF_background1850,3);

BGmask=mean(nSIF_background1850,3);
landmask=BGmask<355;

yesInds=find(landmask);

numRun=numel(yesInds);
total=1;
incr=100;

itter=1;

fl='#!/bin/bash \n';
sl='cd /home/barnhark/seaIceEmergence \n';
tl1='matlab -r "seaIceEmergeDLM_v1(';
tl2='), exit" -nodesktop -nosplash \n';

while total<numRun
    
    start=total;
    stop=min(total+incr-1, numRun);
    
    total=total+incr;
    
    tl=[tl1 num2str(start) ','  num2str(stop) tl2];
   %num2str(itter,'%0.3d')
    
    fileID = fopen(['runDLM_' num2str(itter) '.sh'],'w');
    fprintf(fileID,fl);
    fprintf(fileID,sl);
    fprintf(fileID,tl);
    fclose(fileID);
    
    itter=itter+1;
    
end