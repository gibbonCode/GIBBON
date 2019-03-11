function [varargout]=ellipsoidFit_centered(varargin)

% function [M,s,R,MU]=ellipsoidFit_centered(X,MU)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/11/11
%------------------------------------------------------------------------

%%

switch nargin
    case 1
        X=varargin{1};
        MU=mean(X,1); %Point set mean
    case 2
        X=varargin{1};
        MU=varargin{2};    
end

if isempty(MU)
    MU=mean(X,1); %Point set mean
end

%%

%Centre on mean
X=X-MU(ones(size(X,1),1),:); %Centre points around mean

%Compute singular value decomposition to get principal directions
[~,D,V]=svd(X,0);

%Compute input radii
r = sqrt(sum(X.^2,2));

%Create sphere coordinates
Xr=(V\X')'; %Rotate back
Xc=(D\Xr')'; %Scale back to unit sphere
[theta,phi,~] = cart2sph(Xc(:,1),Xc(:,2),Xc(:,3));
[Xc(:,1),Xc(:,2),Xc(:,3)]=sph2cart(theta,phi,ones(size(theta)));

%Transforming sphere to ellipsoid according to fit
X_fit=Xc;
X_fit=(D*X_fit')'; %Scale
X_fit=(V*X_fit')'; %Rotate
[~,~,r_fit] = cart2sph(X_fit(:,1),X_fit(:,2),X_fit(:,3));

w=mean(r./r_fit); %Determine SVD weight factor

D=D*w; %Scale singular values using weight to convert to ellipsoid semi-axis constants

%Identity
I=eye(4,4);

%Rotation
R=I;
R(1:3,1:3)=V;

%Stretch factors

s=diag(D)';

S=I;
S(1:3,1:3)=D; %Stretch


T=I;
T(1:3,4)=MU(:); %Translation

M  = T * R * S; %The transformation matrix

%%
varargout{1}=M;
varargout{2}=s;
varargout{3}=R;
varargout{4}=MU;

%%


 
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
