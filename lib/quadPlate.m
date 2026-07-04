function [F,V]=quadPlate(varargin)

% function [F,V]=quadPlate(plateDim,plateEl)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2020/07/27
%------------------------------------------------------------------------

%% Parse input 

switch nargin
    case 1
        boxDim=varargin{1};
        boxEl=[10 10 10];        
    case 2
        boxDim=varargin{1};
        boxEl=varargin{2};                
end

%%

dX=boxDim(1); 
dY=boxDim(2); 
nX=boxEl(1); 
nY=boxEl(2);

[X,Y] = meshgrid(linspace(-dX/2,dX/2,nX+1),linspace(-dY/2,dY/2,nY+1)); %Grid

[F,V] = grid2patch(X,Y,zeros(size(X))); %Convert to patch data (quadrilateral faces)

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
