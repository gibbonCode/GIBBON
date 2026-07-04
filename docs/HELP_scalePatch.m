%% scalePatch
% Below is a demonstration of the features of the |scalePatch| function

%% Syntax
% |[Fs,Vs]=scalePatch(F,V,scalefactor);|

%% Description
% This function seperates the nodes (if shared) for all faces. If a
% constant or spatially varying scalefactor is provided the faces are
% scaled (around their mean) as well. 

%% Examples

%%
clear; close all; clc;

%%
% Plot settings
fontSize=15;
plotColor1=0.25.*ones(1,3);
plotColor2=0.75.*ones(1,3);
edgeWidth=2;
markerSize=25;

%% Example: Seperate and scale faces homogeneously
%Defining geodesic dome triangulation
r=1; %sphere radius
n=2; %Refinements
[F,V,~]=geoSphere(n,r);

%Detach nodes and shrink faces
scaleFactor=0.5;
[Fs,Vs]=scalePatch(F,V,scaleFactor);

%%
%Plotting results

cFigure;
subplot(1,2,1); hold on;
gpatch(F,V,'rw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(2);

subplot(1,2,2); hold on;
gpatch(F,V,'rw','none',0.5);
gpatch(Fs,Vs,'bw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(2);

drawnow;

%% Example: Seperate and shrink faces in a spatially varying way

%%
% Create a spatially varying shrink factor between 0 and 1
scaleFactor=V(:,2);
scaleFactor=mean(scaleFactor(F),2);
scaleFactor=scaleFactor-min(scaleFactor(:));
scaleFactor=scaleFactor./max(scaleFactor(:));

%%
% Detach nodes and shrink faces
[Fs,Vs]=scalePatch(F,V,scaleFactor);

%%
% Plotting results

cFigure;
subplot(1,2,1); hold on;
gpatch(F,V,'rw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(2);

subplot(1,2,2); hold on;
gpatch(F,V,'rw','none',0.5);
gpatch(Fs,Vs,'bw','k',1,edgeWidth);
axisGeom(gca,fontSize);
camlight headlight;
view(2);

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
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
