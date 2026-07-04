function [varargout]=diamondLattice(varargin)

% function [Ep,Et,VT,Ct]=diamondLattice(sampleSize,nRepeat,strutThickness)
% ------------------------------------------------------------------------
%
% 2023/04/28 Updated to be a proper diamond lattice, rather than a
% rhombicdodecahedron lattice, and now uses hex2rdl instead. 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        sampleSize=varargin{1};
        nRepeat=[];
        strutThickness=[]; 
        latticePhaseType=[];
    case 2
        sampleSize=varargin{1};
        nRepeat=varargin{2};
        strutThickness=[];  
        latticePhaseType=[];
    case 3
        sampleSize=varargin{1};
        nRepeat=varargin{2};
        strutThickness=varargin{3};
        latticePhaseType=[];
    case 4
        sampleSize=varargin{1};
        nRepeat=varargin{2};
        strutThickness=varargin{3};
        latticePhaseType=varargin{4};        
end

if isempty(nRepeat)
    nRepeat=3;
end

if isempty(strutThickness)
    strutThickness=(sampleSize/nRepeat)/4;        
end

if isempty(latticePhaseType)
    latticePhaseType=1;     
end

%%
% Create a box with hexahedral elements
sampleDimensions=sampleSize*ones(1,3); %Dimensions
sampleElementNumbers=nRepeat*ones(1,3); %Number of elements in each dir
outputStructType=2; %A structure compatible with mesh view
[meshStruct]=hexMeshBox(sampleDimensions,sampleElementNumbers,outputStructType);

%Access elements, nodes, and faces from the structure
Eh=meshStruct.elements; %The elements 
V=meshStruct.nodes; %The nodes (vertices)

F=element2patch(Eh,V);
d=mean(patchEdgeLengths(F,V));

shrinkFactor=strutThickness./(d.*(sqrt(2)./2));

inputStruct.latticePhaseType=latticePhaseType; % 1 = "bubble" centred, 2 = vertex centred, 3 = nested
inputStruct.latticeType=2; % rhombic-dodecahedron (1) or Diamond (2)
inputStruct.shrinkFactor=shrinkFactor; 
inputStruct.removeNonConnected=1;

[Ep,Et,Vs]=hex2rdl(Eh,V,inputStruct); %Create the mesh 

varargout{1}=Ep;
varargout{2}=Et;
varargout{3}=Vs;
if nargout==4
    varargout{4}=[];
    warning('Fourth output no longer supported');
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
