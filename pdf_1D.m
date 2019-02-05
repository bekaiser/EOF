function [ mu sigma minimum maximum var skew flat ] = pdf_1D( u,plotflag,nbins,threshold)
% Calculates statistics for a 1D or 2D array excluding NaNs
% Bryan Kaiser
% 10/18/15

% This function calculates the statistical properties 
% of a given set of random variables "u". "plotflag"
% indicates whether to plot or not (1=yes, 0=no).
% "nbins" is then number of histogram bins, and 
% "threshold" is a positive/negative value which 
% numbers greater/less than will be neglected 
% and converted to NaNs.

% u can be vector or N dimensional matrix

% Set "threshold" to infinity if necessary, it was 
% created for oceanographic data sets with bogus
% values for grid points on land.

% plotflag = 0 => no plot
% plotflag = 1 => pdf plot, no comparison with Gaussian distribution
% plotflag = 2 => pdf plot, comparison with Gaussian distribution

%============================

% Set values above/below a threshold to NaN 
% and find mins/maxes
u(u<-threshold) = NaN;
u(u>threshold) = NaN;

% Remove NaN values (above/below threshold or
% from in the data set:
flagnum = abs(isnan(u)-1);
locnum = find(flagnum);
u = u(locnum); % u is now a vector.

% Compute std dev, skewness, flatness factors:
N = length(u); 
nsigma = zeros(1,N); nskew = nsigma; nflat = nsigma;
mu = mean(u); 
for i = 1:N
       nsigma(i) = (u(i)-mu)^2;
       nskew(i) = (u(i)-mu)^3;
       nflat(i) = (u(i)-mu)^4;
end
var = (sum(nsigma))/N;  % variance
sigma = sqrt(var); % standard deviation
skew = (sum(nskew)/N)/(sigma^3); % skewness
flat = (sum(nflat)/N)/(sigma^4)-3; % flatness
minimum = min(u); % min
maximum = max(u); % max

if plotflag >= 1 % Histogram plot

 % Data distribution   
[counts,centers] =hist(u,nbins);
counts_norm = counts./N;
% counts = number in each bin
% centers = bin center values
% counts_norm = % of total number of samples

if plotflag == 2
% Gaussian distribution
f = zeros(1,nbins); % f = pdf of the normal distribution
for i = 1:nbins
    f(i) = (exp((-(centers(i)-mu).^2)./(2*var))).*(1/(sigma*sqrt(2*pi)));
end
f = f/sum(f);
end

% Plot
figure;
%bar(centers,counts_norm);hold on
if plotflag == 2
plot(centers,f,'r'); % 
end
%xlabel('$$x$$','interpreter','latex','FontSize',20) 
%ylabel('$$P(x)$$','interpreter','latex','FontSize',20) 
%grid on

if plotflag == 2
% Gaussianity check:
Check_pdf_sums_to_one = sum(counts_norm) % pdf sums to 1
Check_Gpdf_sums_to_one = sum(f) % pdf sums to 1

end

end

