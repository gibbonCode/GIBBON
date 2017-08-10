%% subCurve
% Below is a demonstration of the features of the |subCurve| function

%% Syntax
% |[VN]=subCurve(Vt,np,closeLoopOpt);|

%% Description
% The |subCurve| function can be used to increase the point density of the
% input curve by adding evenly spaced points between each of the curve
% segments. 

%% Examples

clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=15;
markerSize1=45;
lineWidth1=2;
lineWidth2=5;
lineWidth3=2;
faceAlpha=0.5;

%% Example: Linearly upsample an open ended curve with intermediate points

%%
% Simulating a curve
Vt=[0 0 0; 10 0 0; 5 10 0; 10 0 10; 0 10 10; ];

%% 
% Upsampling the curve
np=3; %Number of desired intermediate points to be added
[VN]=subCurve(Vt,np); %Using subcurve to upsample curve

%%
% Plotting results
hf1=cFigure;
title('A linearly upsampled curve','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Vt,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(VN,'r.-','lineWidth',lineWidth1/2,'MarkerSize',markerSize1/2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%% Example: Linearly upsample a closed curve with intermediate points

%% 
% Upsampling the curve
np=3; %Number of desired intermediate points to be added
closeLoopOpt=1; %Enable closed loop option
[VN]=subCurve(Vt,np,closeLoopOpt); %Using subcurve to upsample curve

%%
% Plotting results
hf1=cFigure;
title('A linearly upsampled curve with closed end condition','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

plotV(Vt,'k.-.','lineWidth',lineWidth1,'MarkerSize',markerSize1);
plotV(VN,'r.-','lineWidth',lineWidth1/2,'MarkerSize',markerSize1/2);

axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;


%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
 
%% <-- GIBBON footer text --> 
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
