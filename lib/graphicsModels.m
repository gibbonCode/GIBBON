function [F,V]=graphicsModels(varargin)

% function [F,V]=graphicsModels(modelID)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/04/25 Added to GIBBON
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        modelID=1;
    case 1
        modelID=varargin{1};
end

%%
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');

switch modelID
    case {1,'stanford_bunny'}
        fileName=fullfile(pathName,'stanford_bunny_closed.mat');
        meshData=load(fileName);
    case {2,'utah_teapot'}
        fileName=fullfile(pathName,'utah_teapot.mat');
        meshData=load(fileName);
    case {3,'cow'}
        fileName=fullfile(pathName,'cow.mat');
        meshData=load(fileName);
    case {4,'parasaurolophus'}
        fileName=fullfile(pathName,'parasaurolophus.mat');
        meshData=load(fileName);
    case {5,'femur'}
        fileName=fullfile(pathName,'femur.mat');
        meshData=load(fileName);
    case {6,'hip_implant'}
        fileName=fullfile(pathName,'hip_implant.mat');
        meshData=load(fileName);
    case {7,'elephant'} 
        fileName=fullfile(pathName,'elephant.mat');
        meshData=load(fileName);
    case {8,'dolphin'}
        fileName=fullfile(pathName,'dolphin.mat');      
        meshData=load(fileName);
    case {9,'david'}
        % Reduced density version of a David model from SMK - Statens Museum for Kunst
        fileName=fullfile(pathName,'david.mat');
        meshData=load(fileName);
end

F=meshData.F;
V=meshData.V;
 
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
