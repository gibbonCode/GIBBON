function [meshOutput]=importTETGEN(loadNameStruct)

%%

dispStartTitleGibbonCode('Importing TetGen files');

%% IMPORT ELE file
try
    [~,E,elementMaterialID]=importEleFile_tetGen(loadNameStruct.loadName_ele);
catch
%     disp([loadNameStruct.loadName_ele,' import unsuccesful']);
    E=[];
    elementMaterialID=[];
end

%% IMPORT NODE file
try
    [~,V]=importNodeFile_tetGen(loadNameStruct.loadName_node);
    
%     %Check voronoi nodes
%     [pathName,fileName,extension]=fileparts(loadNameStruct.loadName_node);
%     fileName=fullfile(pathName,[fileName,'.v',extension]);
%     if exist(fileName,'file')==2
%         [~,V_vor]=importNodeFile_tetGen(fileName);
%     else
%         V_vor=[];        
%     end
catch
    warning([loadNameStruct.loadName_node,' import unsuccesful']);
    V=[];
    V_vor=[];        
end

%% IMPORT FACES file

try
    [~,F,faceBoundaryID]=importFaceFile_tetGen(loadNameStruct.loadName_face);    
    
%     %Check voronoi faces
%     [pathName,fileName,extension]=fileparts(loadNameStruct.loadName_face);
%     fileName=fullfile(pathName,[fileName,'.v',extension]);
%     if exist(fileName,'file')==2
%         [~,F_vor,faceBoundaryID_vor]=importFaceFile_tetGen(fileName);
%     else
%         F_vor=[];
%         faceBoundaryID_vor=[];
%     end
catch
%     warning([loadNameStruct.loadName_face,' import unsuccesful']);
    F=[];
    faceBoundaryID=[];
    F_vor=[];
    faceBoundaryID_vor=[];
end

%% CONVERT ELEMENTS TO FACES
if ~isempty(E)
    [FE,faceMaterialID]=element2patch(E,elementMaterialID,'tet4');
else
    FE=[];
    faceMaterialID=[];
end

%% Create meshOutput structure
meshOutput.nodes=V; 
meshOutput.facesBoundary=F; 
meshOutput.boundaryMarker=faceBoundaryID;
meshOutput.faces=FE; 
meshOutput.elements=E; 
meshOutput.elementMaterialID=elementMaterialID; 
meshOutput.faceMaterialID=faceMaterialID; 
meshOutput.loadNameStruct=loadNameStruct; 

% meshOutput.facesVoronoi=F_vor; 
% meshOutput.nodesVoronoi=V_vor; 
% meshOutput.boundaryMarkeVoronoi=faceBoundaryID_vor; 

%%
disp(['--- Done --- ',datestr(now)]);
 
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
