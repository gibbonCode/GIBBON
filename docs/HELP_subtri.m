%% subtri
% Below is a demonstration of the features of the |subtri| function

%% Syntax
% |[Fs,Vs]=subtri(F,V,n,uniqueOpt);|

%% Description
% The |subtri| function enables refinement of triangulated data

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

%% Refining a triangle

V=[0 0 0; 1 0 0; 0.5 sqrt(3)/2 0];
F=[1 2 3];

n=0:1:3; %Number of added edge nodes
pColors=gjet(numel(n));

cFigure; 
for q=1:1:numel(n)
    [Fs,Vs]=subtri(F,V,n(q)); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' added edge nodes'],'FontSize',fontSize);
    gpatch(Fs,Vs,pColors(q,:),'k');
    plotV(Vs,'k.','markerSize',markerSize2); 
    plotV(V,'k.','markerSize',markerSize);
    
    axis equal; axis tight; view(2);    
end
drawnow; 

%% Refining a tetrahedron

[V,F]=platonic_solid(1,1);

n=0:1:3; %Number of added edge nodes
pColors=gjet(numel(n));

cFigure; 
for q=1:1:numel(n)
    [Fs,Vs]=subtri(F,V,n(q)); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' added edge nodes'],'FontSize',fontSize);
    gpatch(Fs,Vs,pColors(q,:),'k');
    plotV(Vs,'k.','markerSize',markerSize2); 
    plotV(V,'k.','markerSize',markerSize);    
    axisGeom(gca,fontSize);
    camlight headlight;    
end
drawnow; 

%% Refining triangulated surfaces in general

[F,V]=geoSphere(1,1);

n=[0 1 2 3]; %Number of added edge nodes
pColors=gjet(numel(n));

cFigure; 
for q=1:1:numel(n)
    [Fs,Vs]=subtri(F,V,n(q)); 
    subplot(2,2,q); hold on;
    title([num2str(n(q)),' added edge nodes'],'FontSize',fontSize);
    gpatch(Fs,Vs,pColors(q,:),'k');    
    axisGeom(gca,fontSize);
    camlight headlight;     
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
