function [meshStruct]=tetMeshBox(boxDim,pointSpacing)


%% Create surface mesh
[F,V,faceBoundaryMarker]=triBox(boxDim,pointSpacing);

%%  Mesh model using tetrahedral elements using tetGen

[regionA]=tetVolMeanEst(F,V); %Volume for regular tets

toleranceLevel=1*10^(numOrder(pointSpacing)-8); %tetgen's default is 1e-8
toleranceLevelString=sprintf('%g',toleranceLevel);
stringOpt=['-pq1.2AaY','T',toleranceLevelString];

inputStruct.stringOpt=stringOpt;
inputStruct.Faces=F;
inputStruct.Nodes=V;
inputStruct.holePoints=[];
inputStruct.faceBoundaryMarker=faceBoundaryMarker; %Face boundary markers
inputStruct.regionPoints=mean(V,1); %region points
inputStruct.regionA=regionA;
inputStruct.minRegionMarker=2; %Minimum region marker


[meshStruct]=runTetGen(inputStruct); %Run tetGen 


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
