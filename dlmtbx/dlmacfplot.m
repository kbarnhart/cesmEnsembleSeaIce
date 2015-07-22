function out = dlmacfplot(dlm,ind)
% plot acf of dlm output

if nargin<2, ind = 1; end
p = size(dlm.F,1); % number of series

a = acfnan(dlm.resid2(:,ind));
x = 0:length(a)-1;

colo = [0 0.5 1];

h = plot(x,a,'o-','color',colo,'markerfacecolor',colo);
xlim([x(1),x(end)]);
yl = ylim;ylim([yl(1),1]);
grid
xlabel('lag')
ylabel('acf')
title('Estimated autocorrelation function of the DLM residuals')

% approximate confidence limit, 95%
hh = hline(2/sqrt(dlm.nobs(ind)));
set(hh,'linestyle','--','color','black');
hh = hline(-2/sqrt(dlm.nobs(ind)));
set(hh,'linestyle','--','color','black');

if nargout>0
  out=h;
end
