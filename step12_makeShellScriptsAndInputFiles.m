path='/Volumes/Pitcairn/seaicePPF/northernHemisphere/analysisOutput/';
dirList=dir(strcat(path, 'justNSIF_ensembleAndBG*.nc'));
itter=1;

for f=1:length(dirList)
    fin=strcat(path,dirList(f).name);
    
    nSIF_ensemble = ncread(fin,'nSIF_ensemble') ; 
    nSIF_background1850 = ncread(fin,'nSIF_background1850') ; 
    time_ensemble = ncread(fin,'time_ensemble') ; 
    time_background = ncread(fin,'time_background') ; 
    
    splitStr=strsplit(fin, '.');
    
    saveOutFN=strcat('SIdata.', splitStr{2}, '.', splitStr{3}, '.mat');
    
    save(saveOutFN, 'nSIF_ensemble', 'nSIF_background1850','time_ensemble','time_background')
    
    load (saveOutFN)

    NI=size(nSIF_ensemble,1);
    NJ=size(nSIF_ensemble,2);
    NM=size(nSIF_ensemble,4);
    NY=size(nSIF_ensemble,3);
    NBY=size(nSIF_background1850,3);

    BGmask=mean(nSIF_background1850,3);
    landmask=BGmask<365;

    yesInds=find(landmask);

    numRun=numel(yesInds);
    total=1;
    incr=200;

    
    fl='#!/bin/bash \n';
    sl='cd /home/barnhark/seaIceEmergence \n';
    tl1='matlab -r "seaIceEmergeDLM_v1(';
    tl2='), exit" -nodesktop -nosplash \n';

    while total<numRun

        start=total;
        stop=min(total+incr-1, numRun);

        total=total+incr;

        tl=[tl1 num2str(start) ','  num2str(stop) ,', ''' ,saveOutFN, '''', tl2];
       %num2str(itter,'%0.3d')
%        fileID = fopen(['shellScripts/runDLM.' splitStr{2} '.' splitStr{3} '_' num2str(itter) '.sh'],'w');

        fileID = fopen(['shellScripts/runDLM_' num2str(itter) '.sh'],'w');
        fprintf(fileID,fl);
        fprintf(fileID,sl);
        fprintf(fileID,tl);
        fclose(fileID);

        itter=itter+1;

    end

end