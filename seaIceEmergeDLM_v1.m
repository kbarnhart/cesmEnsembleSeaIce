    
function seaIceEmergeDLM_v1(startyi, stopyi, SIdata)

addpath(genpath('mcmcstat'))
addpath(genpath('dlmtbx'))

% fin='/Volumes/Pitcairn/seaicePPF/northernHemisphere/cesmleOutput/justNSIF_ensembleAndBG.nc';
% nSIF_ensemble = ncread(fin,'nSIF_ensemble') ; 
% nSIF_background1850 = ncread(fin,'nSIF_background1850') ; 
% time_ensemble = ncread(fin,'time_ensemble') ; 
% time_background = ncread(fin,'time_background') ; 

load(SIdata)

splitStr=strsplit(SIdata, '.');

NI=size(nSIF_ensemble,1);
NJ=size(nSIF_ensemble,2);
NM=size(nSIF_ensemble,4);
NY=size(nSIF_ensemble,3);
NBY=size(nSIF_background1850,3);

BGmask=mean(nSIF_background1850,3);
landmask=BGmask<360;

yesInds=find(landmask);

runInds=yesInds(startyi:stopyi);
[nis,njs] = ind2sub(size(landmask),runInds);
nsam_run = 200;
nsimu_run = 1000;


%%
numRuns=numel(runInds);

year_out=zeros(numRuns, 5);
levelmean_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)));
levelstd_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)));
slopemean_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)));
slopestd_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)));

tic

for i=1:numRuns
nii = nis(i);
njj = njs(i);

%% test sites. 
% %DP
% nii=208+1;
% njj=63+1;

% % center of ocean
% nii=193+1;
% njj=77+1;

% Siberian Coast
% nii=156+1;
% njj=60+1;

% nii=133+1;
% njj=70+1;

% % Admunsun
% nii=225+1;
% njj=68+1;

% % iceland sie
% nii=32+1
% njj=81+1


% %kamcha
% nii=174+1;
% njj=24+1;

if landmask(nii, njj)==1
    try 
        tic
        time_ensemble_run=time_ensemble;
        nsam=nsam_run;
        nsimu=nsimu_run;
        fprintf(['nii:',num2str(nii), ' njj: ', num2str(njj) ,'(', num2str(i), '/', num2str(numRuns), ')', '\n'])
        BGdata=reshape(nSIF_background1850(nii,njj,:), [1,NBY]);
        data=reshape(nSIF_ensemble(nii, njj, :,:), [NY, NM]);
        [shiftYear, emergeYear2, emergeYear, levelmean, levelstd, slopemean, slopestd]=runDLM(data(72:end,:), BGdata, time_ensemble_run(72:end), nsam, nsimu);
        toc

        year_out(i,:)=[nii njj shiftYear emergeYear2 emergeYear];
        levelmean_out(i,:)=[nii njj levelmean'];
        levelstd_out(i,:)=[nii njj levelstd'];
        slopemean_out(i,:)=[nii njj slopemean'];
        slopestd_out(i,:)=[nii njj slopestd'];

        %periodically save the files.
        if rem(i, 10)==1
            csvwrite(['output_v3/year_',num2str(startyi),'.csv'], year_out)
            csvwrite(['output_v3/levelmean_',num2str(startyi),'.csv'], levelmean_out)
            csvwrite(['output_v3/slopemean_',num2str(startyi),'.csv'], slopemean_out)
            csvwrite(['output_v3/levelstd_',num2str(startyi),'.csv'], levelstd_out)
            csvwrite(['output_v3/slopestd_',num2str(startyi),'.csv'], slopestd_out)
        end
    catch
        disp 'failed'
    end
end

end
%writeCSV
% csvwrite(['/home/barnhark/seaIceEmergence/year_',num2str(startyi),'.csv'], year_out)
% csvwrite(['/home/barnhark/seaIceEmergence/levelmean_',num2str(startyi),'.csv'], levelmean_out)
% csvwrite(['/home/barnhark/seaIceEmergence/slopemean_',num2str(startyi),'.csv'], slopemean_out)
% csvwrite(['/home/barnhark/seaIceEmergence/levelstd_',num2str(startyi),'.csv'], levelstd_out)
% csvwrite(['/home/barnhark/seaIceEmergence/slopestd_',num2str(startyi),'.csv'], slopestd_out)

%writeCSV
outfile='output_revisions_withbothRCP_nhsh';

csvwrite([outfile,'/year.', splitStr{2}, '.', splitStr{3},'.', num2str(startyi),'.csv'], year_out)
csvwrite([outfile,'/levelmean.', splitStr{2}, '.', splitStr{3},'.', num2str(startyi),'.csv'], levelmean_out)
csvwrite([outfile,'/slopemean.', splitStr{2}, '.', splitStr{3},'.', num2str(startyi),'.csv'], slopemean_out)
csvwrite([outfile,'/levelstd.', splitStr{2}, '.', splitStr{3},'.', num2str(startyi),'.csv'], levelstd_out)
csvwrite([outfile,'/slopestd.', splitStr{2}, '.', splitStr{3},'.', num2str(startyi),'.csv'], slopestd_out)



end
