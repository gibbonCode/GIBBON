%% vectorOrthogonalPair
% Below is a demonstration of the features of the |vectorOrthogonalPair| function

%%
clear; close all; clc;

%% Syntax
% |[a,d]=vectorOrthogonalPair(f);|

%% Description 
% Based on the input vector f this function generates the normalized output
% vectors a and d which are orthogonal to eachother and to f. 
% This function can be useful for describing local element axis systems
% based e.g. in FEBio (whereby the input vector defaults to e3). 

%% Examples

%% Example 1: Creating a triplet of mutually orthogonal vectors 

%%
% Creating example vectors

P=eye(3,3); %Vector origins
V=euler2DCM([0.25*pi 0.25*pi 0.25*pi]); %Input Vectors, rotated directions
% V=eye(3,3); %Input Vectors, x,y,z axes

%%
% Compute mutually orthogonal sets using |vectorOrthogonalPair|
[a,d]=vectorOrthogonalPair(V);

%%
% Visualize the sets

cFigure; 
title('Visualized vectors sets, input vectors=transparent, orthogonal=non-transparent');
quiverVec(P,V,2,'b','none',1,0.2);
quiverVec(P,a,1,'r','none',1,1);
quiverVec(P,d,1,'g','none',1,1);
% quiverVec(P,d,1,'k','none',1,1);
axisGeom; 
camlight headlight; 
drawnow; 

%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
