function [result] = matfilt(fs, fcut, order, data, lowORhigh)

%fs = sampling frequency; fcut= cutoff frequency of filter; order = filter order; lowORhigh = type of filter. 

%switch lowORhigh
%    case 'low'
%        type='low'
%    case 'high'
%        type='high'
%end
type=lowORhigh;
npasses = 3;                       %bi-directional filter
C = (2^(1/npasses) - 1)^0.25;      %the correction factor for number of passes required
Wn1=fcut /(fs/2)/C;%Jay set up
Wn=tan(pi*fcut/fs)/C;%coefficients for a Butterworth or a critically
% damped filtermatlabs; normalaises cutoff frequency


[b,a] = butter(order,Wn,type);    
     
b = double(b);
a = double(a);

[n,m] = size(data);

for i=1:m
  result(:,i) = filtfilt(b,a,data(:,i));
end
