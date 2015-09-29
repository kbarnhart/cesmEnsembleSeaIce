function [shiftYear,emergeYear2, emergeYear, levelmean, levelstd, slopemean, slopestd]=runDLM_logitTest(data, BGdata, time_ensemble, nsam, nsimu)
% written by K. Barnhart
% modified from example code given by Marco Laine associated with the DLM matlab
% toolbox.

% this version of runDLM was specifically modifed to test
% using a logit transformation. 


%% LOGIT TRANSFORM

BGsave=BGdata;
dataSave=data;

BGmeanUT=nanmean(BGdata);
BGstdUT=nanstd(BGdata);
BGubUT=BGmeanUT+2*BGstdUT;
BGlbUT=BGmeanUT-2*BGstdUT;

BGmean=log(max(1,BGmeanUT)/365)-log(1-(max(1,BGmeanUT)/365));
BGstd=log(max(1,BGstdUT)/365)-log(1-(max(1,BGstdUT)/365));
BGub=log(max(1,BGubUT)/365)-log(1-(max(1,BGubUT)/365));
BGlb=log(max(1,BGlbUT)/365)-log(1-(max(1,BGlbUT)/365));

BGdata=log(max(1,BGsave)/365)-log(1-(max(1,BGsave)/365));
data=log(max(dataSave,1)/365)-log(1-(max(dataSave,1)/365));

maxVal=log(364/365)-log(1-(364/365))+2;

%%
% calculate the mean and upper and lower bounds of the background control runDLM

% this is a holdover from a time when I only used a selection of the entire time frame
dataSel=data;
timeSel=time_ensemble;

%%
time = timeSel; % time in years from 1850 2100


% create vectors for the mean and standard deviation of the number of open water days through time.
y = nanmean(dataSel,2); % mean of number of sea ice free days is the
s = nanstd(dataSel, 0,2);

s(isnan(s))=BGstd; %on the strange chance that the std is a nan, use the bg values

% step 1, construct a dlm fit for the interannual variance.
% this will allow us to have a smooth uncertaintyfor the main DLM.
% KRB note: I tried using the standard devation as the uncertainty in the main
% DLM and it produced highly variable prediction intervals. Thus we chose to use
% a smoothed version of the standard deviation.
optionsVar = struct('order',1, 'mcmc',1,'winds',[0 1]); %set DLM options
optionsVar.nsimu=nsimu;
optionsVar.p=1;

% option "winds" specifies which unique parameters should be estimated. here we ask
% set winds=[0 1 ] so that the model uses MCMC to estimate a unique value for
% each of sigma_nu and sigma_beta
sscale=nanmean(s);
sss=s/sscale; % scale for stability
sss_uncert=0.1*mean(sss)*ones(size(sss)); % set uncertainty at 10%

% set prior esimate values for the uncertainty terms sigma_nu and sigma_beta
% (diagonal of parameter W) these values are just initial estimate, because we
% also specify the option "mcmc" and winds =[0 1] which allows for estimation of
% unique values for each of these terms.
w0var = [0 0.1];  % set prior uncertainty for sigma_nu to 0 and sigma_beta to 10%

% inital guess for values of level and trend at time 1
x0var = [mean(sss) 0]; % set uncertanty for sigma_nu to the scaled std and trend
%trend to zero

% construct DLM fit
dlmVar = dlmfit(sss,sss_uncert,w0var,x0var,[],[],optionsVar);

% select the level value and rescale. This will now be used as the uncertainty in
% the main DLM model
news=dlmVar.x(1,:)'*sscale;
newss=news;


%% Step 2 construct main DLM fit.
% scale y, (typically the number of sea ice free days per year) for stability
% by dividing by the standard deviation.
ys = stdnan(reshape(y, 1, numel(y)));
yy = y./ys;
% also scale the standard deviation for stability
ss = newss./ys;
ss(ss<0.01)=0.01; % variance of zero  -> nonstable, set at 1 percent.

%%
% Prior means for some components of |W|, the model error matrix.
wtrend = 0.02; % set prior estimate for the uncertanty in trend (sigma_beta)
% KRB note: investigated values from 0.02 to 0.5 - didn't have a big effect.
w0 = [0 wtrend];
% again, this is just an initial guess - because we set the term "mcmc" and winds=[0 1]
% within the model we estimate unique values for sigma_nu and sigma_beta using mcmc

% set inital guesses for the level and trend at time 1. use the normalized mean of the background
% for the level and zero for the trend.
x0 = [BGmean/ys 0];% 0]; % initial values, use standardize background mean values

%%
% Calculate the DLM smoother solution, do MCMC over some components in the matrix |W|.

% for the sea ice component, our DLM has two components
% a mean state
% a linear trend (specified by stating "order"=1)

% at test sites autocovariance structure looks good.
% KRB note: tested including an AR(1) term, which had little effect.

options = struct('order',1, 'mcmc',1,'nsimu',1000,'winds',[0 1]);
options.p=size(y,2);
dlm = dlmfit(yy,ss,w0,x0,[],[],options);

%%
% Produce sample from the model states using |dlmsmosam|. It accounts the posterior uncertainty
% in W using the MCMC chain in |dlm.chain|.
 % number of sampled to draw from the posterior
dlm_sample = dlmsmosam(dlm,nsam);
% %%

%%
% Sample trend statistics form DLM sample
levelsamp = squeeze(dlm_sample(1,:,:)); % sample of levels
levelmean = mean(levelsamp, 2);              % their mean
levelstd = std(levelsamp, 0, 2);

% exploration of a 10 year-long (+-5) running trend.
off=10;
t2 = mean((levelsamp(off:end,:)-levelsamp(1:end-off+1,:))/off,2); % mean trend
s2 = std((levelsamp(off:end,:)-levelsamp(1:end-off+1,:))/off, 0,2);

% figure(8); clf
% hold on
% plot(timeSel(off:end),t2+2*s2)
% plot(timeSel(off:end),t2)
% plot(timeSel(off:end),t2-2*s2, 'r');grid;

slopesamp = squeeze(dlm_sample(2,:,:));
slopemean=mean(slopesamp,2);
slopestd=std(slopesamp, 0, 2);


%% find years (current implementation allows for not shifting and not emerging)
% neg=(t2+(2*s2)<0) ;
% pos=(t2-(2*s2)>0) ;

inside=(t2+(2*s2)>0).*(t2-(2*s2)<0);
shiftInds=find(inside==0)+off;

% choose first of longest time.
split1=find(diff(shiftInds)>1);
if ~isempty(split1)
    split=[1; split1; numel(shiftInds)];
    for k = 1:numel(split)-1
        sizes(k)=split(k+1)-split(k);
    end
    bigger=find(sizes==max(sizes), 1, 'last');
    startsplit=split(bigger+1);
    stopslit=split(bigger+1);
    biggerInds=shiftInds(startsplit:stopslit);
    shiftInd=min(biggerInds);

else
   shiftInd=min(shiftInds);
end

% if numel(t2)==max(shiftInds)
%    % assume that the system has reached full open water conditons.
%
%    dumInd=find(diff(shiftInds)>1, 1, 'last');
%    shiftInds(dumInd:end)=[];
% end
% shiftInd=max(shiftInds)+off+1; % choose last one, + one because we are choosing the last year inside.


if isempty(shiftInd)
    shiftYear=1800;
elseif shiftInd==1
    shiftYear=1800;
else
    shiftYear=time(shiftInd);
end

% some of this needs to be different in logit-land
% 1) first, use the mean and std of bg BEFORE transform to construct CIs
% then transform
% 2) transform totally out threhsold. 


over2=(levelmean+(2*levelstd)<BGlb);
under2=(levelmean-(2*levelstd)>BGub);

totalOutThreshUT=363;
totalOutThresh=log(max(1,totalOutThreshUT)/365)-log(1-(max(1,totalOutThreshUT)/365));

totallyOut=levelmean>totalOutThresh;

out2=over2+under2+totallyOut;
emergeInd2=find(out2>0, 1, 'first');
if isempty(emergeInd2)
    emergeYear2=2110;
else
    emergeYear2=time(emergeInd2);
end


over=(levelmean+(2*news)<BGlb);
under=(levelmean-(2*news)>BGub);

% there seem to be some problems with the sea ice present all the time
% areas,


out=over+under+totallyOut;

emergeInd=find((out>0)&(levelmean>3), 1, 'first');
if isempty(emergeInd)
    emergeYear=2110;
else
    emergeYear=time(emergeInd);
end


%% plots
% figure;
% dlmplotfit(dlm, time, ys)
% title('Title');xlabel('time');ylabel('Open Water Days Per Year]')
%
% figure;
% dlmplotdiag(dlm, time, ys)

% redo plot QQ plot for split up residuals
figure;
dlmplotdiag_logit_krbMOD(dlm, time, ys)


% %%%%%
% figure;
% hold on
% for i=1:5:nsam
%   plot(time,squeeze(dlm_sample(1,:,i)),'-')
% end
% hold off
% 
% figure;
% hold on
% for i=1:5:nsam
%   plot(time,squeeze(dlm_sample(2,:,i)),'-')
% end
% hold off
%%%%%%

% % The next figure shows prior and posterior distributions for
% % standard deviations from the diagonal of model error matrix |W|.
% figure(5); clf
% mcmcplot(dlm.chain,[],dlm.res,'denspanel',2);
% subplot(2,1,1);title('prior and posterior for variance parameters');xlabel('parameter w(2,2)')
% subplot(2,1,2);title('');xlabel('parameter w(3,3)')

figure; clf
confband(timeSel(off:end),t2,s2);grid;
hold on
xlim([time(1),time(end)]); % match axis to other plots
title('Yearly Trend');
ylabel('Days Per Year')
hold off

figure; clf
plot(timeSel, dataSel, 'k')%, 'Alpha',0.2)
hold on
confband(timeSel,levelmean,levelstd);grid;
% get at prediction interval
plot(timeSel, levelmean+2*news ,'r--')
plot(timeSel, levelmean-2*news,'r--')

plot([1920, 2120],[BGmean, BGmean])
plot([1920, 2120],[BGub, BGub])
plot([1920, 2120],[BGlb, BGlb])

plot([shiftYear, shiftYear],[0,10])
plot([emergeYear, emergeYear],[0,10])
plot([emergeYear2, emergeYear2],[0,10])

hold off

xlim([time(1),time(end)+20]); % match axis to other plots
title('Model Output, Logit Transformation');
ylabel('Logit Transformation: Days Per Year')
xlabel('Time')


% %% UNLOGIT Transform (dlm_sample)
% 
% dlm_sample_save=dlm_sample;
% 
% % to return to days of open water, first scale by ys, then un-logit
% % transform, then multiply the level term by 365
% 
% dlm_sample=(exp(ys.*dlm_sample_save))./(1+exp(ys.*dlm_sample_save));
% dlm_sample(1,:,:)=365.*dlm_sample(1,:,:);
% 
% 
% BGmeanUT
% BGubUT
% BGlbUT

end
