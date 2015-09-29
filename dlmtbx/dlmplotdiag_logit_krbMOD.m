function out = dlmplotdiag_logit_krbMOD(dlm, t, ys, yind)
%DLMPLOTDIAG  Default plot for DLM model diagnostics
% dlmplotdiag(dlm, t, ys, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% ys - optional y axis scale factor for plots (not used)
% yind - which series to plot, default = 1

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $
% revised by Katy Barnhart 2015/09
% to show subsetted residuals 

if nargin<2
  t = (1:size(dlm.y,1))'; % yscale for plots
end
if nargin<3
  ys = 1; % yscale for plots
end
if nargin<4
  yind = 1; % which obs column to plot
end

% select where 
bottomThresh=3;
topThresh=362;

bottomThreshLogit=log(bottomThresh/365)-log(1-(bottomThresh/365));
topThreshLogit=log(topThresh/365)-log(1-(topThresh/365));

subInd1=find((ys.*dlm.y)<=bottomThreshLogit);
subInd2=find(((ys.*dlm.y)>bottomThreshLogit).*((ys.*dlm.y)<topThreshLogit));
subInd3=find((ys.*dlm.y)>=topThreshLogit);

% mutliple qqplots
subplot(3,1,1)
dlmqqplot_KRBmod(dlm,yind, subInd1);
title('Normal probability plot for subset residuals (nSIF<=3)')
subplot(3,1,2)
dlmqqplot_KRBmod(dlm,yind,subInd2);
title('Normal probability plot for subset residuals (nSIF>3 and nSIF<362)')
subplot(3,1,3)
dlmqqplot_KRBmod(dlm,yind, subInd3);
title('Normal probability plot for subset residuals (nSIF>=362)')

if nargout>1
  out = [];
end
