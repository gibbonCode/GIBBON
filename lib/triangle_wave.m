function y=triangle_wave(x)




p=2.*pi;
y=(4.*(abs((x./p)-floor((x./p)+0.5))))-1;

