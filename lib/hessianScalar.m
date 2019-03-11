function [H]=hessianScalar(varargin)

% function [H]=hessianScalar(U,v,cellOpt)
% ------------------------------------------------------------------------
%
% This function computes the Hessian matrix of the scalar function U for
% each point in U. U may be a vector, a 2D matrix or a 3D matrix. The
% vector v denotes the points spacing between the data entries in U. If v
% is not supplied the spacing is assumed to be homogeneous and unity. If
% the input is n dimensional array consisting of m entries then the output
% is a matrix (if cellOpt==0) the size of mx(n^2) (whereby the colum
% entries define the entries in a Hessian matrix and row entries relate to
% elements in the input array). If cellOpt==1 then the output is reformed
% into a cell array that matches the size of the input aray. Each cell
% entry then contains the nxn Hessian matrix. 
%
% See also: |gradient|,|hessian|,|jacobian|,|cellEig|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/10 Updated
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        U=varargin{1};
        v=[];
        cellOpt=0;
    case 2
        U=varargin{1};
        v=varargin{2};
        cellOpt=0;
    case 3
        U=varargin{1};
        v=varargin{2};
        cellOpt=varargin{3};
end
%%

nDims=ndims(U);
if isvector(U)
    nDims=1;
end

if isempty(v)
    v=ones(1,nDims); %Use default
end

if nDims<=3    
    switch nDims
        case 1
            %Compute first order derivative
            [dUdx] = gradient(U,v);
            
            %Compute second order derivatives
            [dUdxdx] = gradient(dUdx,v);
            
            H_mat=dUdxdx(:);
        case 2
            %Compute first order derivative
            [dUdx,dUdy] = gradient(U,v(1),v(2));
            
            %Compute second order derivatives
            [dUdxdx,dUdxdy] = gradient(dUdx,v(1),v(2));
            [dUdydx,dUdydy] = gradient(dUdy,v(1),v(2));
            
            H_mat=[dUdxdx(:) dUdxdy(:)...
                dUdydx(:) dUdydy(:)];                                 
        case 3            
            %Compute first order derivative
            [dUdx,dUdy,dUdz] = gradient(U,v(1),v(2),v(3));
            
            %Compute second order derivatives
            [dUdxdx,dUdxdy,dUdxdz] = gradient(dUdx,v(1),v(2),v(3));
            [dUdydx,dUdydy,dUdydz] = gradient(dUdy,v(1),v(2),v(3));
            [dUdzdx,dUdzdy,dUdzdz] = gradient(dUdz,v(1),v(2),v(3));
            
            H_mat=[dUdxdx(:) dUdxdy(:) dUdxdz(:)...
                   dUdydx(:) dUdydy(:) dUdydz(:)...
                   dUdzdx(:) dUdzdy(:) dUdzdz(:)];
    end
    
    if cellOpt==1
        H_mat=reshape(H_mat',nDims,nDims,size(H_mat,1));
        
        %Cell array of Hessian matrices
        H=reshape(mat2cell(H_mat,nDims,nDims,ones(size(H_mat,3),1)),size(U));
    else
        H=H_mat;
    end
else     
    warning('Function only defined for 1D, 2D or 3D scalar valued arrays');
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
