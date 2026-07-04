%% export_STL_txt
% Below is a demonstration of the features of the |export_STL_txt| function

%%
clear; close all; clc;

%%
% Plot settings
faceAlpha=0.5;
fontSize=25; 

%% Get patch data sets
[F1,V1]=stanford_bunny;
meanV=mean(V1,1);
V1=V1-meanV(ones(1,size(V1,1)),:);
V1=V1./max(V1(:))/2;
[F2,V2]=parasaurolophus;
meanV=mean(V2,1);
V2=V2-meanV(ones(1,size(V2,1)),:);
V2=V2./max(V2(:));
V2=V2+max(V1,[],1);

%% Create the stlStruct

stlStruct.solidNames={'stanford_bunny','parasaurolophus'}; %names of parts
stlStruct.solidVertices={V1,V2}; %Vertices
stlStruct.solidFaces={F1,F2}; %Faces
stlStruct.solidNormals={[],[]}; %Face normals (optional)

%%
% Plotting the models 
pColors=gjet(numel(stlStruct.solidNames));

cFigure;
title('Patch data to export to multi-solid STL','fontSize',fontSize);

for q=1:1:numel(stlStruct.solidNames)
    F=stlStruct.solidFaces{q}; %Faces
    V=stlStruct.solidVertices{q}; %Vertices
    
    gpatch(F,V,pColors(q,:));
end

axisGeom
camlight('headlight');
lighting phong; axis off; 
drawnow;

%% Exporting an STL file from the multi-solid patch data

%Set main folder and fileName
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','STL'); 
fileName=fullfile(pathName,'stanford_bunny_multi.stl'); 

export_STL_txt(fileName,stlStruct);

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
