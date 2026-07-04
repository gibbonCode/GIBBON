%% hex2tet
% Below is a demonstration of the features of the |hex2tet| function

%% Syntax
% |[TET,Vtet,C]=hex2tet(HEX,V,C,tetOpt);|

%% Description
%

%%
clear; close all; clc;

%%
% Plot settings
fontSize=25;
faceAlpha1=0.25;
edgeWidth=3;
markerSize=75;
cMap=gjet(6);

%% Examples
%

%% Converting a hexahedral element to tetrahedral elements

%%
% Creating an example hexahedral element
V=[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1;]; %nodes
E=1:8; %Element

%%

[F,CF]=element2patch(E);  %Patch data for plotting

cFigure;
subplot(2,3,1); hold on;
title('Original element set','FontSize',fontSize);
gpatch(F,V,cMap(1,:),'k',1,edgeWidth);
% patchNormPlot(F,V,0.25);
plotV(V,'k.','MarkerSize',markerSize);
colormap(cMap);
axisGeom(gca,fontSize);
axis off; 

for q=1:1:5
    
    % Subdeviding the hexahedral element
    convertMethod=q; %Corresponse
    [Es,Vs,Cs]=hex2tet(E,V,[],convertMethod);
    [Fs,CFs]=element2patch(Es); %Patch data for plotting
    
    subplot(2,3,q+1); hold on;
    title(['Converted, method: ',num2str(q)],'FontSize',fontSize);
    gpatch(Fs,Vs,cMap(q+1,:),'k',0.5,edgeWidth);
    patchNormPlot(Fs,Vs,0.25);
    plotV(Vs,'k.','MarkerSize',markerSize);
    
    colormap(cMap);
    axisGeom(gca,fontSize);
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
