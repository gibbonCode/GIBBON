function [Vg]=sweepCurveBezier(varargin)

% function [Vg]=sweepCurveBezier(p1,p2,n1,n2,numSteps,f)
% ------------------------------------------------------------------------
%
%
%
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 4
        p1=varargin{1};
        p2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numPoints=10;
        f=1/3;
    case 5
        p1=varargin{1};
        p2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numPoints=varargin{5};    
        f=1/3;        
    case 6    
        p1=varargin{1};
        p2=varargin{2};
        n1=varargin{3};
        n2=varargin{4};
        numPoints=varargin{5};
        f=varargin{6};            
end

n1=vecnormalize(n1);
n2=vecnormalize(n2);

d=sqrt(sum(p1-p2).^2); %Distance between input points
w=d*f; %"tangency" weights
Pb=[p1; p1+w*n1;  p2-w*n2; p2]; %Bezier points

Vg=bezierCurve(Pb,numPoints);

Vg=evenlySampleCurve(Vg,numPoints,'pchip',0);
 
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
