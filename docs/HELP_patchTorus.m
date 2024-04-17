%% patchTorus
% Below is a demonstration of the features of the |patchTorus| function

%%
clear; close all; clc;

%% Syntax
% |[F,V]=patchTorus(r,nr,rc,nc,patchType);|

%% Description
% This function generates patch for a torus with desired dimensions and
% mesh type

%% Examples

%%
% Plot settings
fontSize=10;
faceColor='b';
faceAlpha=1;
edgeColor='k';
edgeWidth=2;
markerSize=25;

%% Building a mesh of a torus

%Torus parameters
r=1; %Sphere radius
rc=2.5; %Central radius
nr=16;
nc=25;
patchTypes={'quad','tri_q','tri','honey'};

% Open figure for plotting
cFigure;

%Plot the various mesh types
pColors=gjet(4);
for q=1:1:4
    [F,V]=patchTorus(r,nr,rc,nc,patchTypes{q});
    subplot(2,2,q); hold on;
    title([patchTypes{q}],'FontSize',fontSize,'Interpreter','none');
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hp=patch('Faces',F,'Vertices',V);
%     [hn]=patchNormPlot(F,V,0.3);
    set(hp,'FaceColor',pColors(q,:),'FaceAlpha',faceAlpha,'lineWidth',edgeWidth,'edgeColor',edgeColor);
    camlight headlight;
    set(gca,'FontSize',fontSize);
    view(3); axis tight;  axis equal;  grid on;
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
