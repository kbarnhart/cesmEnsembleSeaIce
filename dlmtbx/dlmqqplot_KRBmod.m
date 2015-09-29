function out=dlmqqplot_KRBmod(dlm,ind, ind2)
%DLMQQPLOT normal probability plot for DLM fit residuals 

if nargin<2, ind = 1; end
% p = size(dlm.F,1); % number of series

if ~isempty(ind2)

    y = dlm.resid2(ind2,ind);
    igood = not(isnan(y));
    yy = y(igood);
    h=qqplot(yy);
    %set(h,'markerfacecolor','black');
    grid on;

    if nargout>0
      out=h;
  
end
end
