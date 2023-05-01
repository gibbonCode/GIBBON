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
