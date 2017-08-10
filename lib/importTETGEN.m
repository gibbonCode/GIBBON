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
catch
    warning([loadNameStruct.loadName_node,' import unsuccesful']);
    V=[];
end

%% IMPORT FACES file
try
    [~,F,faceBoundaryID]=importFaceFile_tetGen(loadNameStruct.loadName_face);
catch
%     warning([loadNameStruct.loadName_face,' import unsuccesful']);
    F=[];
    faceBoundaryID=[];
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

%%
disp(['--- Done --- ',datestr(now)]);

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
