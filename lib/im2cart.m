function [X,Y,Z]=im2cart(I,J,K,v)

% function [X Y Z]=im2cart(I,J,K,v)
% ------------------------------------------------------------------------
% This function converts the image coordinates I,J,K to the cartesian
% coordinates X,Y,Z using the voxel dimension v. 
%
% I,J,K can be scalars, vectors or matrices. 
% v is a vector of length 3 where v(1), v(2) and v(3) correspond to the
% voxel dimensions in the x,y and z direction respectively. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% ------------------------------------------------------------------------

X=(J-0.5).*v(2);
Y=(I-0.5).*v(1);
Z=(K-0.5).*v(3);