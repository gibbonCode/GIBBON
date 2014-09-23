function [H]=hessianScalar(U,v)

% function [H]=hessianScalar(U,v)
% ------------------------------------------------------------------------
%
% This function computes the Hessian matrix of the scalar function U for
% each point in U. U may be a vector, a 2D matrix or a 3D matrix. The
% vector v denotes the points spacing between the data entries in U. If v
% is not supplied the spacing is assumed to be homogeneous and unity. For
% 1D data the output is a vector with dimensions size(U). If the input is
% 2D or 3D the output is a cell array with dimensions size(U) containing
% 2x2 or 3x3 matrices respectively. 
%
% See also: |gradient|,|hessian|,|jacobian|,|cellEig|
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/11/2013
%------------------------------------------------------------------------


nDims=ndims(U);

if nDims<=3
    
    if nargin==1
        v=ones(1,nDims);
    end
    
    if isvector(U)
        nDims=1;
    end
    
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
    
    H=H_mat;
%     H_mat=reshape(H_mat',nDims,nDims,size(H_mat,1));
%     
%     %Cell array of Hessian matrices
%     H=reshape(mat2cell(H_mat,nDims,nDims,ones(size(H_mat,3),1)),size(U));
else     
    warning('Function only defined for 1D, 2D or 3D scalar valued arrays');
end




