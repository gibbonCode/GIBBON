%% quad2tri
% Below is a demonstration of the features of the |quad2tri| function

%% Syntax
% |[varargout]=quad2tri(varargin);|

%% Description 
% 
%% Examples 
% 

%%
clear; close all; clc;

%%
% Plot settings
fontSize=10;

%% Example: Converting a quandrangulated surface to a triangulated surface

%%
% Create example geometry

[X,Y,Z]=peaks(15);
[F,V,CV]=surf2patch(X,Y,Z,Z);
[C]=vertexToFaceMeasure(F,CV);

%%
% Visualisation

cFigure;
title('Quadrangulation','FontSize',fontSize);
gpatch(F,V,C,'k');
axisGeom(gca,fontSize);
colormap(gjet(250)); colorbar;
camlight('headlight'); 
drawnow; 

%%
% Convert triangular faces to quadrilateral faces

cFigure;
subplot(2,3,1);
title('Quadrangulation','FontSize',fontSize);
gpatch(F,V,C,'k');
axisGeom;
view(-30,70);
camlight headlight
colormap(gjet(250)); 

optionSet={'f','b','x','e','a'};
titleString={'forw. slash','backw. slash','cross','edge length','angle'};
numOptions=numel(optionSet);
for q=1:1:numOptions
    [Ft,Vt,Ct]=quad2tri(F,V,optionSet{q},C);
    subplot(2,3,q+1);
    title(['Option ''',optionSet{q},''' ',titleString{q}],'FontSize',fontSize);
    gpatch(Ft,Vt,Ct,'k');
    axisGeom;
    view(-30,70);
    camlight headlight
    colormap(gjet(250));
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
