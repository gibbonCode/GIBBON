function [varargout]=patchMarchDistMapIterative(varargin)


%%

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        indStart=1;
        W=ones(size(V,1),1);
        optStruct.waitBarOn=0;
        optStruct.numIterationsMax=200;
        optStruct.toleranceLevel=1e-3;
    case 3
        F=varargin{1};
        V=varargin{2};
        indStart=varargin{3};
        W=ones(size(V,1),1);
        optStruct.waitBarOn=0;
        optStruct.numIterationsMax=200;
        optStruct.toleranceLevel=1e-3;
    case 4
        F=varargin{1};
        V=varargin{2};
        indStart=varargin{3};
        W=varargin{4};
        optStruct.waitBarOn=0;
        optStruct.numIterationsMax=200;
        optStruct.toleranceLevel=1e-3;
    case 5
        F=varargin{1};
        V=varargin{2};
        indStart=varargin{3};
        W=varargin{4};
        optStruct=varargin{5};
end

if isempty(W)
    W=ones(size(V,1),1);
end

if isempty(indStart)
    indStart=1;
end

if isfield(optStruct,'numIterationsMax')
    numIterationsMax=optStruct.numIterationsMax;
else
    numIterationsMax=200;
end

if isfield(optStruct,'toleranceLevel')
    toleranceLevel=optStruct.toleranceLevel;
else
    toleranceLevel=1e-3;
end

%%

i = F(:);
j = [F(:,2); F(:,3); F(:,1)];
k = [F(:,3); F(:,1); F(:,2)];
numPoints = size(V,1);
x  = V(i,:);
x1 = V(j,:) - x;
x2 = V(k,:) - x;

Inv1 = @(M,d)[M(:,4)./d -M(:,3)./d -M(:,2)./d M(:,1)./d];
Inv  = @(M)Inv1(M, M(:,1).*M(:,4) - M(:,3).*M(:,2));
Mult = @(M,u)[M(:,1).*u(:,1) + M(:,3).*u(:,2)  M(:,2).*u(:,1) + M(:,4).*u(:,2)];

U = getoptions(optStruct, 'U', []);
if isempty(U)
    U = zeros(numPoints,1);
end

% inner product matrix
C = [dot(x1,x1,2) dot(x1,x2,2) dot(x2,x1,2) dot(x2,x2,2)];

S = Inv(C);

% a = <S*1,1>
a = sum(S,2);

w = W(i);

% edge length
L1 = sqrt(dot(x1,x1,2));
L1 = L1(:).*w(:);
L2 = sqrt(dot(x2,x2,2));
L2 = L2(:).*w(:);

E=nan(numPoints,1);

ii=1:1:numel(i);

for q=1:numIterationsMax
    uj = U(j);
    uk = U(k);
    u = [uj uk];
    b = sum([S(:,1)+S(:,3)  S(:,2)+S(:,4)].*u,2);
    c = sum(Mult(S,u).*u,2)-w.^2;
    delta = max( b.^2 - a.*c, 0);
    d = (b + sqrt(delta) )./a; % solution

    alpha = Mult(S,u-repmat(d, 1,2)); % direction of the update

    logic_d =alpha(:,1)>0 | alpha(:,2)>0;

    % update along edges
    d1 = L1 + uj(:);
    d2 = L2 + uk(:);
    d = d(:);
    d(logic_d) = min(d1(logic_d), d2(logic_d));

    %Create U1 using sparse array to get minima
    maxd=max(d(:));
    df=abs(maxd-d);
    U1 = sparse(i,ii,df,numPoints,max(ii),numel(d));
    U1 =maxd-full(max(U1,[],2));
    U1(U1==0) = Inf;
    U1(indStart) = 0; % boundary condition

    % enforce monotony
    if min(U1-U)<-1e-5
      warning('Monotony problem');
    end

    % Store error
    if q==1
        e1=norm(U-U1,'fro');
    end

    E(q) = norm(U-U1,'fro')/e1;
    if E(q)<toleranceLevel
        break;
    end

    % Update
    U = U1;
end

%Output
varargout{1}=U;
varargout{2}=E;
 
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
