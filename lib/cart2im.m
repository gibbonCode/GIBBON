function [varargout]=cart2im(varargin)

% function [I J K]=cart2im(X,Y,Z,v)
% ------------------------------------------------------------------------
% This function converts the cartesian coordinates X,Y,Z to image
% coordinates I,J,K using the voxel dimension v.
%
% X,Y,Z can be scalars, vectors or matrices. 
% v is a vector of length 3 where v(1), v(2) and v(3) correspond to the
% voxel dimensions in the x,y and z direction respectively. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% 2019/05/24 Added handling of variable input and output
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        XYZ=varargin{1};         
        v=[];        
        compactInput=1;
    case 2
        XYZ=varargin{1};        
        v=varargin{2};        
        compactInput=1;
    case 3
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        v=[];
        compactInput=0;
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        v=varargin{4};
        compactInput=0;
end

%Get X,Y,Z if provided in a single array
if compactInput==1
    X=XYZ(:,1);
    Y=XYZ(:,2);
    if size(XYZ,2)>2
        Z=XYZ(:,3);
    else
        Z=[];
    end
end

if isempty(v)
    v=ones(1,3);
end

if numel(v)==1
    v=v*ones(1,3);
end

if isempty(Z)
    Z=0.5.*v(3)*ones(size(X)); %Initialize as "mid-slice" for 2D input
end

%% Process conversion
I=(Y./v(1))+0.5;
J=(X./v(2))+0.5;
K=(Z./v(3))+0.5;

%% Collect output
switch nargout
    case 1
        varargout{1}=[I(:) J(:) K(:)];
    otherwise
        varargout{1}=I;
        varargout{2}=J;
        varargout{3}=K;
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
