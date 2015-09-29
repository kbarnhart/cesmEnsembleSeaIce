function out = dlmplotdiag_krbMOD(dlm, t, ys, yind)
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

subInd1=find((ys.*dlm.y)<=bottomThresh);
subInd2=find(((ys.*dlm.y)>bottomThresh).*((ys.*dlm.y)<topThresh));
subInd3=find((ys.*dlm.y)>=topThresh);

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
