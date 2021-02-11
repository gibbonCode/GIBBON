function [N]=randUnitVec(n)

% function [N]=randUnitVec(n)
%-------------------------------------------------------------------------
% This function generates points on a unit sphere
%
% 2020/12/22 Created and added to GIBBON
%-------------------------------------------------------------------------

%%

phi=rand_angle([n,1]); % range [0 2*pi>, similar to 2*pi*rand(n,1);
theta=acos((2*rand(n,1))-1); 
N=[cos(phi).*sin(theta) sin(phi).*sin(theta) cos(theta)]; %Unit vectors

