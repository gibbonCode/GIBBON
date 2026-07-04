function [F,V]=import_obj_geom(fileName)

% function [F,V]=import_obj_geom(fileName)
% -----------------------------------------------------------------------
% This function imports only the geometry data from the obj file defined by
% fileName. A single body is assumed. 
%
% Change log: 
% 2021/04/30 KMM: Fixed bug relating to / symbols for faces 
% 2021/04/30 KMM: Speeded up by using single cellfun based loop
% 2023/06/01 KMM: Improved speed, added handling of mixed face types 
% -----------------------------------------------------------------------

%% Import text into string array

S_OBJ = readlines(fileName);

%% Get logic for different text components

% Find lines to exclude
Lc=contains(S_OBJ,'# ');
L_mtl=contains(S_OBJ,'mtllib');
L_mat=contains(S_OBJ,'usemtl');
L_exc = Lc | L_mtl | L_mat; %Exclude logic

%Create vertex and face logics
Lv=contains(S_OBJ,'v ') & ~L_exc;
Lf=contains(S_OBJ,'f ') & ~L_exc;

%% Get vertices

nVertices = nnz(Lv);
V=zeros(nVertices,3);
s=S_OBJ(Lv);
for q=1:1:nVertices
    V(q,:)=sscanf(s(q),'v %f %f %f')';
end

%% Get faces

nFaces = nnz(Lf); 
FC=repmat({zeros(1,3)},nFaces,1);  %Allocate cell array for triangles (with vertex/normal data) for the moment (multiple face types may occur)
nf=zeros(nFaces,1);
s=S_OBJ(Lf);
sf=regexprep(s,'/+[0-9]+|f', ' '); %Remove f and get just face vertex indices
for q=1:1:nFaces
    FC{q}=sscanf(sf(q),'%d')';    
    nf(q)=numel(FC{q});
end

faceTypeSet=unique(nf);
numFaceTypes=numel(faceTypeSet);
F=cell(numFaceTypes,1);
for q=1:1:numFaceTypes
    logicNow = nf==faceTypeSet(q);
    F{q}=cell2mat(FC(logicNow));
end

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
