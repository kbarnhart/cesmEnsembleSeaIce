% run dlm for monthly extents
close all
clear all

addpath(genpath('mcmcstat'))
addpath(genpath('dlmtbx'))


load /Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/OWExtent_analysis.mat 

time=1921:2100;
nsam = 200;
nsimu = 1000;
areaConversion=1000000*1000*1000;
tic
pK={'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sept', 'Oct', 'Nov', 'Dec'};

for i=1:12
    
    BGdata=BGmonthMean(72:end,i)'./areaConversion;
    data=squeeze(iceMonthMean(:,72:end,i))'./areaConversion;
    
    [shiftYear(i), emergeYear2(i), emergeYear(i), levelmean, levelstd, slopemean, slopestd]=...
    runDLM_extent(data, BGdata, time, nsam, nsimu, pK{i});

    toc

end