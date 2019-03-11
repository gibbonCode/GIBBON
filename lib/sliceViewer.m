function hf=sliceViewer(varargin)


%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
        viewerType=1;
    case 2
        M=varargin{1};
        v=varargin{2};
        viewerType=1;
    case 3
        M=varargin{1};
        v=varargin{2};
        viewerType=varargin{3};
end
        
%% Start viewer

switch viewerType
    case 1
        hf=sv2(M,v);
    case 2
        hf=sv3(M,v);
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
