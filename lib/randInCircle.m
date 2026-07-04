function [V]=randInCircle(varargin)

% function [V]=randInCircle(n,R)
% ------------------------------------------------------------------------
% Create random and uniform distribution of n points in a circle with
% radius R. The output array V is an nx2 2D coordinate array. 
%
% 2021/06/28 Created KMM
% ------------------------------------------------------------------------

%%

switch nargin
    case 1
        n=varargin{1};
        R=1;
    case 2
        n=varargin{1};
        R=varargin{2};
    otherwise
        error('Wrong number of input arguments. Define 1 or 2 input arguments'); 
end

%%

r = R.*sqrt(rand(n,1)); %Radii
theta = rand_angle([n,1]); %rand(n,1) * 2*pi; %Angles
V=[r.*cos(theta) r.*sin(theta)]; %Points

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
