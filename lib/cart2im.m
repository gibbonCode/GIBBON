function [I,J,K]=cart2im(X,Y,Z,v)

% function [I J K]=cart2im(X,Y,Z,v)
% ------------------------------------------------------------------------
% This function converts the cartesian coordinates X,Y,Z to image
% coordinates I,J,K using the voxel dimension v.
%
% X,Y,Z can be scalars, vectors or matrices. 
% v is a vector of length 3 where v(1), v(2) and v(3) correspond to the
% voxel dimensions in the x,y and z direction respectively. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% ------------------------------------------------------------------------

I=(Y./v(1))+0.5;
J=(X./v(2))+0.5;
K=(Z./v(3))+0.5;