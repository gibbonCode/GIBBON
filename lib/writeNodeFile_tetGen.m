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

cell2txtfile(nodeFileName,T,0,0);
dispDoneGibbonCode;
 
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
