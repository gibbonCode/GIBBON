function [varargout]=isInsideTR(varargin)



%%

switch nargin
    case 1
        TR=varargin{1};
        QP=TR.Points; 
        ti=[];
        toleranceMagnitude=0;
    case 2
        TR=varargin{1};
        QP=varargin{2};
        ti=[];
        toleranceMagnitude=0;
    case 3
        TR=varargin{1};
        QP=varargin{2};
        ti=varargin{3};
        toleranceMagnitude=0;
    case 4
        TR=varargin{1};
        QP=varargin{2};
        ti=varargin{3};
        toleranceMagnitude=varargin{4};
    otherwise 
        error('Wrong number of input arguments');
end

if isempty(ti)
    ti=1:size(TR.ConnectivityList,1);
end
    
%Get barycentric coordinates of points
baryCoords=cartesianToBarycentric(TR,ti(:),QP);

logicFoundEnclosing=all(baryCoords>=-toleranceMagnitude,2);

varargout{1}=logicFoundEnclosing;
if nargout==2
    varargout{2}=baryCoords;
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
