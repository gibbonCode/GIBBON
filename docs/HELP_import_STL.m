%% import_STL
% Below is a demonstration of the features of the |import_STL| function

%% Syntax
% |[stlStruct] = import_STL(fileName);|

%% Description
% Use |import_STL| to import binary or txt type STL files.

%% Examples

clear; close all; clc;

%%
% Plot settings
fontSize=25; 

%% Import a txt type STL file as patch data

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'femur.stl'); 
[stlStruct] = import_STL(fileName);

F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};

%%
% Merging nodes example
[F,V]=mergeVertices(F,V);

%%
% Plotting the model

cFigure;
title('Imported patch data from STL','fontSize',fontSize);
gpatch(F,V,'g');
axisGeom;
camlight('headlight');
lighting phong; axis off; 
drawnow;

%% Import a binary type STL file as patch data

%Set main folder
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'femur_bin.stl'); 
[stlStruct] = import_STL(fileName);

F=stlStruct.solidFaces{1};
V=stlStruct.solidVertices{1};

%%
% Merging nodes example
[~,ind1,ind2]=unique(pround(V,5),'rows');
V=V(ind1,:);
F=ind2(F);

%%
% Plotting the model

cFigure;
title('Imported patch data from STL','fontSize',fontSize);
gpatch(F,V,'g');
axisGeom;
camlight('headlight');
lighting phong; axis off; 
drawnow;

%%

[M,ellipStretchFit,R_fit,MU]=ellipsoidFit_centered(V,[]);

%Create sphere
[F_fit,V_fit,~]=geoSphere(4,1);

%Transforming sphere to ellipsoid
V_fit_t=V_fit;
V_fit_t(:,end+1)=1;
V_fit_t=(M*V_fit_t')'; %Rotate polyhedron
V_fit=V_fit_t(:,1:end-1);
cFigure; hold on; 
title('The true (green) and fitted ellipsoid (red) and axis directions (solid, transparant respectively)','FontSize',fontSize);

gpatch(F,V,'gw','k',1);
gpatch(F_fit,V_fit,'rw','none',0.2);

axisGeom;
camlight('headlight');
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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
