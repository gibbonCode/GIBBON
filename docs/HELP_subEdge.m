%% subEdge
% Below is a demonstration of the features of the |subEdge| function

%% Syntax
% |[Fs,Vs]=subEdge(F,V,n,uniqueOpt);|

%% Description
% The |subEdge| function enables refinement of the edges of patch

%% Examples

clear; close all; clc;

%% 
% Plot Settings
fontSize=15;
faceAlpha=1;
edgeColor=0.2*ones(1,3);
edgeWidth=1.5;
markerSize=35; 
markerSize2=20; 

%% Refining the edges of a triangle

V=[0 0 0; 1 0 0; 0.5 sqrt(3)/2 0];
F=[1 2 3];

n=0:1:3; %Number of added edge nodes
pColors=gjet(numel(n));
cFigure; 
for q=1:1:numel(n)
    [Fs,Vs]=subEdge(F,V,n(q)); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' added edge nodes'],'FontSize',fontSize);
    gpatch(Fs,Vs,pColors(q,:),edgeColor,faceAlpha);
    plotV(Vs,'k.','markerSize',markerSize2);
    plotV(V,'k.','markerSize',markerSize);    
    set(gca,'FontSize',fontSize);
    view(2); axis tight;  axis equal;  axis off; 
end

%% Refining  the edges of 3D patch data

[F,V]=quadSphere(2,1);

n=0:1:3; %Number of added edge nodes
pColors=gjet(numel(n));
cFigure; 
for q=1:1:numel(n)
    [Fs,Vs]=subEdge(F,V,n(q)); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' added edge nodes'],'FontSize',fontSize);
    gpatch(Fs,Vs,pColors(q,:),edgeColor,faceAlpha); 
    plotV(Vs,'k.','markerSize',markerSize2);
    plotV(V,'k.','markerSize',markerSize);    
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  axis off; 
    camlight headlight; 
end

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
