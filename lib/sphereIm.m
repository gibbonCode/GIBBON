function M=sphereIm(varargin)

% function M=sphereIm(n,logicFlag,nd)
%-------------------------------------------------------------------------
% This function computes a so called "sphere image" where the radius n is
% given 
%
% 
%-------------------------------------------------------------------------

%% Parse input
switch nargin
    case 1
        n=varargin{1};
        logicFlag=1;
        nd=3;
    case 2
        n=varargin{1};
        logicFlag=varargin{2};
        nd=3;
    case 3
        n=varargin{1};
        logicFlag=varargin{2};
        nd=varargin{3};
end

%% Compute image coordinates
nf=floor(n); 

%Set ouput image size
if nd>1 %2D or larger
    siz=((2*nf)+1)*ones(1,nd);
else %1D
    siz=[((2*nf)+1) 1]; %Add singular 2nd dimension
end

%Subscript index array for all voxels
P=ind2subn(siz,1:1:prod(siz));

%% Compute distance image

M=reshape(sqrt(sum((P-nf-1).^2,2)),siz);

%Convert to logic if needed
if logicFlag==1
    M=M<=n;
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
