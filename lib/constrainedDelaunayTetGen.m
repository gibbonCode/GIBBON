function [varargout]=constrainedDelaunayTetGen(V,C)


%% Create TetGen input structure

if ~isempty(C)
    inputStruct.Faces=C; %Add face constraints if not empty
    inputStruct.stringOpt='-pQY';
else
    inputStruct.stringOpt='-Q'; 
end
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=ones(size(C,1),1); %Face boundary markers
inputStruct.regionPoints=[]; %region points

%% Run TetGen
[meshOutput]=runTetGen(inputStruct); %Run tetGen

TR = triangulation(meshOutput.elements,meshOutput.nodes);

switch nargout
    case 1
        varargout{1}=TR;
    case 2
        varargout{1}=TR;
        varargout{2}=meshOutput;
end
 
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
