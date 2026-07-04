function [varargout]=im2cart(varargin)

% function [X Y Z]=im2cart(I,J,K,v)
% ------------------------------------------------------------------------
% This function converts the image coordinates I,J,K to the cartesian
% coordinates X,Y,Z using the voxel dimension v. 
%
% I,J,K can be scalars, vectors or matrices. 
% v is a vector of length 3 where v(1), v(2) and v(3) correspond to the
% voxel dimensions in the x,y and z direction respectively. 
%
% This function maps the row, column, slice coordinates I,J,K to
% "real-world" Cartesian coordinates X,Y,Z based on: 
% X=(J-0.5).*v(2);
% Y=(I-0.5).*v(1);
% Z=(K-0.5).*v(3);
%
% Note that the columns relate to X while rows relate to Y. 
%
% A single coordinate array may also be specified whereby the columns
% define the I, J, and K coordinates. If a single output is requested the 
% output will also consist of such an array whereby columns are the X, Y,
% and Z coordinates
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% 2019/05/24 Added handling of variable output
% 2019/06/25 Added handling of 2D input
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        IJK=varargin{1};         
        v=[];        
        compactInput=1;
    case 2
        IJK=varargin{1};        
        v=varargin{2};        
        compactInput=1;
        I=IJK(:,1);
        J=IJK(:,2);
        if size(IJK,2)>2
            K=IJK(:,3);
        else
            K=[];
        end
    case 3
        I=varargin{1};
        J=varargin{2};
        K=varargin{3};
        v=[];
        compactInput=0;
    case 4
        I=varargin{1};
        J=varargin{2};
        K=varargin{3};
        v=varargin{4};
        compactInput=0;
end

%Get I, J, K if provided in a single array
if compactInput==1
    I=IJK(:,1);
    J=IJK(:,2);
    if size(IJK,2)>2
        K=IJK(:,3);
    else
        K=[];
    end
end

%Initialize K as ones if missing (e.g. for 2D input)
if isempty(K)
   K=ones(size(I)); 
end

%Default voxel size if missing
if isempty(v)
    v=ones(1,3);
end

%Expand voxel size isotropically if single value is provided
if numel(v)==1
    v=v*ones(1,3);
end

%% Process conversion (shift and scale)

X=(J-0.5).*v(2);
Y=(I-0.5).*v(1);
Z=(K-0.5).*v(3);

%% Collect output
switch nargout
    case 1
        varargout{1}=[X(:) Y(:) Z(:)];
    otherwise
        varargout{1}=X;
        varargout{2}=Y;
        varargout{3}=Z;
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
