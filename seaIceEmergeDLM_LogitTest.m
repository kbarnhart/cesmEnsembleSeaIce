  
close all
clear all

addpath(genpath('mcmcstat'))
addpath(genpath('dlmtbx'))

pathOut='assumptionCheckFigures';

SIdata='SIdata.nh.RCP85.mat';
load(SIdata)

splitStr=strsplit(SIdata, '.');

NI=size(nSIF_ensemble,1);
NJ=size(nSIF_ensemble,2);
NM=size(nSIF_ensemble,4);
NY=size(nSIF_ensemble,3);
NBY=size(nSIF_background1850,3);

BGmask=mean(nSIF_background1850,3);
landmask=BGmask<360;

nsam_run = 200;
nsimu_run = 1000;


%%
filename = '/Users/katherinebarnhart/git/cesmEnsembleSeaIce/assumptionCheckPoints.csv';
delimiter = ',';
startRow = 2;

% Format string for each line of text:
%   column1: text (%s)
%	column2: double (%f)
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Allocate imported array to column variable names
sitename = dataArray{:, 1};
nis = dataArray{:, 2};
njs = dataArray{:, 3};


% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

%%
numRuns=numel(nis);

year_out=zeros(numRuns, 5,2);
levelmean_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)),2);
levelstd_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)),2);
slopemean_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)),2);
slopestd_out=zeros(numRuns, 2+numel(time_ensemble(72:end,:)),2);

tic

for i=1:numRuns
nii = nis(i);
njj = njs(i);


if landmask(nii, njj)==1
    try 
        site=sitename{i};
        
        fprintf(site)

        tic
        time_ensemble_run=time_ensemble;
        nsam=nsam_run;
        nsimu=nsimu_run;
        fprintf(['nii:',num2str(nii), ' njj: ', num2str(njj) ,'(', num2str(i), '/', num2str(numRuns), ')', '\n'])
        BGdata=reshape(nSIF_background1850(nii,njj,:), [1,NBY]);
        data=reshape(nSIF_ensemble(nii, njj, :,:), [NY, NM]);
        
        
        close all
        [shiftYearLogit, emergeYear2Logit, emergeYearLogit, levelmeanLogit, levelstdLogit, slopemeanLogit, slopestdLogit]   =   runDLM_logitTest(   data(72:end,:), BGdata, time_ensemble_run(72:end), nsam, nsimu);
        hfigs = get(0, 'children');
        
        for m = 1:length(hfigs)
            filename=[pathOut '/RegLogitTest.' site '.logit.' num2str(m) '.pdf'];                                %Bring Figure to foreground
            if strcmp(filename, '0')                        %Skip figure when user types 0
                continue
            else
                saveas(hfigs(m), [filename '.pdf']) %.pdf
            end
        end
        
        close all
        [shiftYear, emergeYear2, emergeYear, levelmean, levelstd, slopemean, slopestd]                                      =   runDLM(             data(72:end,:), BGdata, time_ensemble_run(72:end), nsam, nsimu);
        hfigs = get(0, 'children');
        
        for m = 1:length(hfigs)
            filename=[pathOut '/RegLogitTest.' site '.untransformed.' num2str(m) '.pdf'];                                %Bring Figure to foreground
            if strcmp(filename, '0')                        %Skip figure when user types 0
                continue
            else
                saveas(hfigs(m), [filename '.pdf']) %.pdf
            end
        end
        toc

        year_out(i,:,1)=[nii njj shiftYear emergeYear2 emergeYear];
        levelmean_out(i,:,1)=[nii njj levelmean'];
        levelstd_out(i,:,1)=[nii njj levelstd'];
        slopemean_out(i,:,1)=[nii njj slopemean'];
        slopestd_out(i,:,1)=[nii njj slopestd'];

        year_out(i,:,2)=[nii njj shiftYearLogit emergeYear2Logit emergeYearLogit];
        levelmean_out(i,:,2)=[nii njj levelmeanLogit'];
        levelstd_out(i,:,2)=[nii njj levelstdLogit'];
        slopemean_out(i,:,2)=[nii njj slopemeanLogit'];
        slopestd_out(i,:,2)=[nii njj slopestdLogit'];
    % save figures for each Location:
    
        
        close all

    
    
    catch
        disp 'failed'
    end
end

end



% compare year chosen:

figure()
hold on
plot(year_out(:,3,1)) % normal, shift
plot(year_out(:,3,2)) % logit transform, shift
plot(year_out(:,5,1)) % normal, emerge
plot(year_out(:,5,2)) % logit transform, emerge
h = gca;
h.XTick =1:numRuns;
h.XTickLabel = sitename;
h.XTickLabelRotation=90;
legend('shift', 'logit shift', 'emerge', 'logit emerge')
hold off

