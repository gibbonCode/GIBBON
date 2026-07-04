function [varargout]=patchBoundary(varargin)

%function [Eb,E,indBoundary]=patchBoundary(F)
%-------------------------------------------------------------------------
%
% 
%
% Change log: 
% 2018/10/08 Expanded to handle cell input for F
% 2021/07/30
% 2021/10/08 Simplified, especially cell handling
% 2021/10/08 No longer needs vertices as input
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        F=varargin{1};
    case 2
        F=varargin{1};
        warning('Second input (vertices) no longer required. Update code to avoid future error.');
end

%% Get boundary edges

% Get non-unique edges
E=patchEdges(F,0);
   
% Get boundary edge indices
Es=sort(E,2); %Sort so edges with same nodes have the same rows
[~,~,~,countUse]=cunique(Es,'rows'); %Count occurances
logicBoundary=countUse==1;

% Boundary edges
Eb=E(logicBoundary,:);

%% Gather output
varargout{1}=Eb; %Boundary edges
varargout{2}=E; %All edges
if nargout==3
    varargout{3}=find(logicBoundary); %Indices for boundary edges
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
