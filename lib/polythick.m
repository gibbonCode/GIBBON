function [Xu Yu Xl Yl]=polythick(x,y,t)

% function [Xu Yu Xl Yl]=polythick(x,y,t)
% ------------------------------------------------------------------------
% This function thickens a polygon defined by x, y. It generates the upper
% (Xu, Yu) and lower (Xl, Yl) coordinates depending on the thickness t. 
%
% The slope is calculated at each coordinate and points are copied upwards
% and downwards (thickening) orthogonal to the local slope.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 27/08/2008
% ------------------------------------------------------------------------

%% 

%Calculating 'derivate' using the function diff
y_diff=diff(y)./diff(x);
x_diff=(x(1:end-1)+x(2:end))./2;

%Interpolating points of derivative at x coordinates of the polygon
y_der = interp1(x_diff,y_diff,x,'cubic');

%Calculating slope angle in radians
alpha_dir=atan(y_der);

%Creating upper and lower coordinates
Xu=x-((0.5*t)*(sin(alpha_dir)));
Yu=y+((0.5*t)*(cos(alpha_dir)));
Xl=x+((0.5*t)*(sin(alpha_dir)));
Yl=y-((0.5*t)*(cos(alpha_dir)));

%% END