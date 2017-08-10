function [varargout]=polyTube(V,optStruct)

d=pathLength(V);

if isfield(optStruct,'r')
    r=optStruct.r;
else
    r=min(diff(d))/2; %Half the minimum point spacing
end

if isfield(optStruct,'C')
    C=optStruct.C;
else
    C=[1:1:size(V,1)]'; %Point indices
end

if isfield(optStruct,'nr')
    nr=optStruct.nr;
else
    nr=3;
end

if isfield(optStruct,'patchType')
    patchType=optStruct.patchType;
else
    patchType='tri';
end

%%

t=linspace(0,2*pi,nr+1);
t=t(1:end-1);
x=r*sin(t);
y=r*cos(t);
z=zeros(size(t));
Vc=[x(:) y(:) z(:)];
V1=V(2:end,:)-V(1:end-1,:);
V1(end+1,:)=V1(end,:);

V1m=(V1(2:end,:)+V1(1:end-1,:))/2;
V1m=[V1(1,:); V1m];

[Fni,Vni,~]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),V1m(:,1),V1m(:,2),V1m(:,3),[],[1 1]);

VC=repmat(Vc,[size(V,1),1]);

vr=[0 0 1]+rand(1,3); vr=vecnormalize(vr);
V3=V1m+vr(ones(size(V1m,1),1),:); V3(:,3)=V3(:,3)+0.5; V3=vecnormalize(V3);
V2=cross(V1m,V3,2); V2=vecnormalize(V2);
V3=cross(V1m,V2,2); V3=vecnormalize(V3);

for q=2:1:size(V1m,1)
    V3(q,:)=V3(q-1,:);
    V2(q,:)=cross(V3(q,:),V1m(q,:)); V2=vecnormalize(V2);
    V3(q,:)=cross(V1m(q,:),V2(q,:)); V3=vecnormalize(V3);
end

X=V1m(:,ones(nr,1))';
Y=V1m(:,2*ones(nr,1))';
Z=V1m(:,3*ones(nr,1))';
V1m_rep=[X(:) Y(:) Z(:)];



X=V2(:,ones(nr,1))';
Y=V2(:,2*ones(nr,1))';
Z=V2(:,3*ones(nr,1))';
V2_rep=[X(:) Y(:) Z(:)];

X=V3(:,ones(nr,1))';
Y=V3(:,2*ones(nr,1))';
Z=V3(:,3*ones(nr,1))';
V3_rep=[X(:) Y(:) Z(:)];

R_mat=[V3_rep V2_rep V1m_rep];
R_mat=R_mat(:,[1 4 7 2 5 8 3 6 9]);
[VC]=vectorTensorProductArray(VC,R_mat);

X=V(:,ones(nr,1))';
Y=V(:,2*ones(nr,1))';
Z=V(:,3*ones(nr,1))';
VC_offset=[X(:) Y(:) Z(:)];

VC=VC+VC_offset;

colorInd=reshape([1:1:size(V1m_rep,1)]',[nr,size(V,1)])';

X=reshape(VC(:,1),[nr,size(V,1)])';
Y=reshape(VC(:,2),[nr,size(V,1)])';
Z=reshape(VC(:,3),[nr,size(V,1)])';

[Fs,Vs,Cs_ind] = surf2patch(X,Y,Z,colorInd);

I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
Fs_sub=sub2ind(size(Z),I,J);
Fs=[Fs;Fs_sub];

switch patchType
    case 'quad'
    case 'tri'
        [Fs,Vs]=quad2tri(Fs,Vs,'f');
end

%Create colors

Cs_d=repmat(d,[1,nr]);
Cs_d=Cs_d(:);
% Cs_d=vertexToFaceMeasure(Fs,Cs_d(:));

Cs_rgb=abs(V1m_rep(Cs_ind,:));

Cs_C=repmat(C,[1,nr]);
Cs_C=Cs_C(:);

%%

switch nargout
    case 2
        varargout{1}=Fs;
        varargout{2}=Vs;
    case 3
        varargout{1}=Fs;
        varargout{2}=Vs;
        varargout{3}=Cs_C;
    case 4
        varargout{1}=Fs;
        varargout{2}=Vs;
        varargout{3}=Cs_C;
        varargout{4}=Cs_rgb;
    case 5
        varargout{1}=Fs;
        varargout{2}=Vs;
        varargout{3}=Cs_C;
        varargout{4}=Cs_rgb;
        varargout{5}=Cs_d;
    otherwise
        error('Wrong number of output arguments');
end
 
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
