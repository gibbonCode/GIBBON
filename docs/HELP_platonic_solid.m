%% platonic_solid
% Below is a demonstration of the features of the |platonic_solid| function

%%
clear; close all; clc;

%%
% Plot settings
fig_color='w'; fig_colordef='white';
fontSize=15;
faceColor='b';
faceAlpha=0.75;
edgeColor='k';
edgeWidth=2;
markerSize=5;

%%

hf=figuremax(fig_color,fig_colordef); % Open figure for plotting
pColor=jet(5);
for q=1:1:5
    %Defining the faces (F) and vertices (V) of a platonic solid
    [V,F]=platonic_solid(q,1); %q indicates solid type, r is the radius
    
    subplot(2,3,q); hold on;
    
    hp=gpatch(F,V,pColor(q,:),'k',faceAlpha);
    
    %Plotting face normals
    [hn]=patchNormPlot(F,V);
    
    axisGeom(gca,fontSize);    
    camlight('headlight'); lighting flat;
    axis off;
end
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
