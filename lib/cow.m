function [F,V]=cow

% [F,V]=cow
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% cow model. 
%
% This model consists of 5804 triangular faces and 2903 vertices.
% 
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/21
%------------------------------------------------------------------------
%%

% filePath=mfilename('fullpath');
% loadPath=fullfile(fileparts(fileparts(filePath)),'data','stl','cow.stl');
% [stlStruct] = import_STL_txt(loadPath);
% F=stlStruct.solidFaces{1};
% V=stlStruct.solidVertices{1};
%  [F,V,ind1,ind2]=mergeVertices(F,V,numDigitsMerge);

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');
fileName=fullfile(pathName,'cow.mat');
D=load(fileName);
F=D.F;
V=D.V;


 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
