function [E]=curveToEdgeList(varargin)

% function [E]=curveToEdgeList(V,closedCurveOpt)
% ------------------------------------------------------------------------
% This function creates an edge array E for the input curve defined by N.
% The input parameter N can be of the following type: 
% * A single scalar
% In this case it is assumed N defines the number of points on the curve
% * An nxm array representing n points and m dimensions, in this case it is
% assumed N represents the vertex array for the curve
% * A row or column array, in this case it is assumed that N defines the
% indices for the points defining the curve. 
% 
% If closeCurveOpt=1 it is assumed the start and end of the curve should be
% attached. 
% 
% 2023/05/31 Updated description and input handling. Added closed curve
% option. 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        N=varargin{1}; 
        closedCurveOpt=0;
    case 2
        N=varargin{1}; 
        closedCurveOpt=varargin{2};
end

%Check input type
if numel(N)==1 %the size of the list is specified
    indList=(1:1:N)';
elseif isvector(N) %The indices are provided
    indList=mcol(N);
else %ordered vertices are provided
   indList=(1:1:size(N,1))';
end

%Deal with closed curve option
if closedCurveOpt
    indList(end+1)=indList(1);
end

%% Create edge array
E=[indList(1:end-1) indList(2:end)];

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
