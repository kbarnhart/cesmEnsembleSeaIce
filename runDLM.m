function [shiftYear,emergeYear2, emergeYear, levelmean, levelstd, slopemean, slopestd]=runDLM(data, BGdata, time_ensemble, nsam, nsimu)
% written by K. Barnhart
% modified from example code given by Marco Laine associated with the DLM matlab
% toolbox.



% calculate the mean and upper and lower bounds of the background control runDLM
BGmean=nanmean(BGdata);
BGstd=nanstd(BGdata);
BGub=BGmean+2*BGstd;
BGlb=BGmean-2*BGstd;

% this is a holdover from a time when I only used a selection of the entire time frame
dataSel=data;
timeSel=time_ensemble;

%%
time = timeSel; % time in years from 1850 2100
y = nanmean(dataSel,2); % mean of number of sea ice free days is the
s = nanstd(dataSel, 0,2);
s(isnan(s))=BGstd; %on the strange chance that the std is a nan, use the bg values

% step 1, construct a dlm fit for the interannual variance.
% this will allow us to get the prediction interval

optionsVar = struct('order',1, 'mcmc',1,'winds',[0 1]); %set DLM options
optionsVar.nsimu=nsimu;
optionsVar.p=1;

sscale=nanmean(s);
sss=s/sscale; % scale for stability
sss_uncert=0.1*mean(sss)*ones(size(sss));

% set values for the
w0var = [0 0.1];  % set uncertainty for sigma_beta to 10%
x0var = [mean(sss) 0]; % set uncertanty for sigma_nu to the scaled std

% construct DLM
dlmVar = dlmfit(sss,sss_uncert,w0var,x0var,[],[],optionsVar);

%
news=dlmVar.x(1,:)'*sscale;
newss=news;

toc
%%
% scale y for stability
ys = stdnan(reshape(y, 1, numel(y)));
yy = y./ys;
ss = newss./ys;
ss(ss<0.01)=0.01; % variance of zero  -> nonstable, set at 1 percent.

%%
% Prior means for some components of |W|, the model error matrix.
wtrend = 0.02; % set prior estimate for the uncertanty in trend (sigma_beta)
% KRB note: investigated values from 0.02 to 0.5 - didn't have a big effect.
w0 = [0 wtrend];
x0 = [BGmean/ys 0];% 0]; % initial values, use standardize background mean values

%%
% Calculate the DLM smoother solution, do MCMC over some components in the matrix |W|.

% for the sea ice component, our DLM has two components
% a mean state
% a linear trend
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
levelsamp = ys*squeeze(dlm_sample(1,:,:)); % sample of levels
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

slopesamp = ys*squeeze(dlm_sample(2,:,:));
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

over2=(levelmean+(2*levelstd)<BGlb);
under2=(levelmean-(2*levelstd)>BGub);
totallyOut=levelmean>363;

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
% figure(1);
% dlmplotfit(dlm, time, ys)
% title('Title');xlabel('time');ylabel('Open Water Days Per Year]')
% %
% figure(2);
% dlmplotdiag(dlm, time, ys)
%
%
% %%%%%
% figure(3);
% hold on
% for i=1:5:nsam
%   plot(time,ys*squeeze(dlm_sample(1,:,i)),'-')
% end
% hold off
%
% figure(4);
% hold on
% for i=1:5:nsam
%   plot(time,ys*squeeze(dlm_sample(2,:,i)),'-')
% end
% hold off
% %%%%%%
%
% % The next figure shows prior and posterior distributions for
% % standard deviations from the diagonal of model error matrix |W|.
% figure(5); clf
% mcmcplot(dlm.chain,[],dlm.res,'denspanel',2);
% subplot(2,1,1);title('prior and posterior for variance parameters');xlabel('parameter w(2,2)')
% subplot(2,1,2);title('');xlabel('parameter w(3,3)')

%
% figure(6); clf
% confband(timeSel(off:end),t2,s2);grid;
% hold on
% xlim([time(1),time(end)]); % match axis to other plots
% title('Yearly Trend');
% ylabel('Days Per Year')
% hold off
%
% figure(7); clf
% plot(timeSel, dataSel, 'k')%, 'Alpha',0.2)
% hold on
% confband(timeSel,levelmean,levelstd);grid;
% % get at prediction interval
% plot(timeSel, levelmean+2*news ,'r--')
% plot(timeSel, levelmean-2*news,'r--')
%
% plot([1920, 2100],[BGmean, BGmean])
% plot([1920, 2100],[BGub, BGub])
% plot([1920, 2100],[BGlb, BGlb])
%
%
% plot([shiftYear, shiftYear],[0,365])
% plot([emergeYear, emergeYear],[0,365])
% plot([emergeYear2, emergeYear2],[0,365])
%
% hold off
%
% xlim([time(1),time(end)]); % match axis to other plots
% title('Yearly Trend');
% ylabel('Days Per Year')

end
