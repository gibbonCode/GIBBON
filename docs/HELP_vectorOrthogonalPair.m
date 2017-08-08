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
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
