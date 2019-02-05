function [ PDF ] = pdf_2D( u,v,plotflag,nbins,threshold )
% Calculates the JPDF for two 1D arrays excluding NaNs
% Bryan Kaiser
% 2/14/16

% Set values above/below a threshold to NaN 
% and find mins/maxes
u(u<-threshold) = NaN;
u(u>threshold) = NaN;
v(v<-threshold) = NaN;
v(v>threshold) = NaN;

x = linspace(min(u),max(u),nbins); 
y = linspace(min(v),max(v),nbins); 
[X,Y] = meshgrid(x,y); 

pdf = hist3([u,v],{x y}); % bivariate histogram
PDF = pdf'./length(u); % prbability density function

if plotflag == 1
    figure;
    surf(X,Y,PDF);
    %set(h,'edgecolor','none');
end

Check_pdf_sums_to_one = sum(sum(PDF)) % pdf sums to 1

end

