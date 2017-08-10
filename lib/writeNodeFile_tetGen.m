function [T]=writeNodeFile_tetGen(inputStruct)

%%

dispStartTitleGibbonCode('Writing NODE file');

%% PARSE INPUT STRUCTURE

V=inputStruct.Nodes;
nodeFileName=inputStruct.modelName; 

%Force extension to be .node
[pathstr,name,~] = fileparts(nodeFileName);
nodeFileName=fullfile(pathstr,[name,'.node']);

%%

V_id=1:1:size(V,1);
V_field=[V_id(:) V];
V_char=sprintf('%i %0.16e %0.16e %0.16e \n',V_field');
V_cell = regexp(V_char, '\n', 'split'); 

if numel(V_cell)>1
    V_cell=V_cell(1:end-1);
end

T={'#PART 1 - Node list';'#num nodes, num dimensions, num attributes, num boundary markers'};
numNodes=size(V,1);
numDims=3; 
numAtr=0;
boundMarker=0;
doubleList=[numNodes numDims numAtr boundMarker];
charList=sprintf('%i %i %i %i',doubleList');
T(end+1)={charList};
T(end+1)={'#Node ID, x, y, z,attribute,boundary marker'};
T(end+1:end+numel(V_cell))=V_cell;

%% SAVING TXT FILE

cell2txtfile(nodeFileName,T,0);
dispDoneGibbonCode;


 
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
