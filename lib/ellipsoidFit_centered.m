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


 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
