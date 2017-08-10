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
