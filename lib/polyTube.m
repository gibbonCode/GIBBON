function [varargout]=polyTube(V,optionStruct)

% function [varargout]=polyTube(V,optionStruct)
% ------------------------------------------------------------------------
%
%
%
%
% Change log (started 2019/03/28): 
% 2019/03/28 Added spatially varying radii
% 2019/03/29 Added capped ends option (for triangulated patch type) 
% ------------------------------------------------------------------------

%% Parse input

numPoints=size(V,1);

%Create default option structure
d=pathLength(V);
defaultOptionStruct.r=min(diff(d))/2; %Half the minimum point spacing
defaultOptionStruct.C=[1:1:numPoints]'; %Point indices
defaultOptionStruct.nr=3;
defaultOptionStruct.patchType='quad';
defaultOptionStruct.closeOpt=0;

%Parse input structure (complete with default)
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1);

%Access options
C=optionStruct.C;
nr=optionStruct.nr;
patchType=optionStruct.patchType;
closeOpt=optionStruct.closeOpt;
r=optionStruct.r;
if numel(r)==1
    r=r.*ones(numPoints,1);
end

%%

t=linspace(0,2*pi,nr+1);
t=t(1:end-1);
x=sin(t);
y=cos(t);
z=zeros(size(t));
Vc=[x(:) y(:) z(:)];
V1=V(2:end,:)-V(1:end-1,:);
V1(end+1,:)=V1(end,:);

V1m=(V1(2:end,:)+V1(1:end-1,:))/2;
V1m=[V1(1,:); V1m];

VC=repmat(Vc,[numPoints,1]);

RC=repmat(r(:)',[size(Vc,1),1]);
RC=RC(:);
VC=VC.*RC(:,ones(1,size(VC,2))); %Scale radii

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

colorInd=reshape((1:1:size(V1m_rep,1))',[nr,numPoints])';

X=reshape(VC(:,1),[nr,numPoints])';
Y=reshape(VC(:,2),[nr,numPoints])';
Z=reshape(VC(:,3),[nr,numPoints])';

[Fs,Vs,Cs_ind] = surf2patch(X,Y,Z,colorInd);

I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))'];
J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
Fs_sub=sub2ind(size(Z),I,J);
Fs=[Fs;Fs_sub];

%% Create colors
Cs_d=repmat(d,[1,nr]);
Cs_d=Cs_d(:);
% Cs_d=vertexToFaceMeasure(Fs,Cs_d(:));

Cs_rgb=vecnormalize(abs(V1m_rep(Cs_ind,:)));

Cs_C=repmat(C,[1,nr]);
Cs_C=Cs_C(:);


%%

switch patchType
    case 'quad'
        if closeOpt
            error('Capping cylinder is currently only supported for the tri patch type');
        end
    case 'tri'
        [Fs,Vs]=quad2tri(Fs,Vs,'f');
        if closeOpt
            indTop=numPoints:numPoints:size(Vs,1);
            [Ft,Vt]=regionTriMesh3D({Vs(indTop,:)},[],0);
            
            indBottom=1:numPoints:size(Vs,1);
            [Fb,Vb]=regionTriMesh3D({Vs(indBottom,:)},[],0);
            
            [Fs,Vs,Cs]=joinElementSets({Fs,Ft,Fb},{Vs,Vt,Vb});
            [Fs,Vs,ind1,ind2]=mergeVertices(Fs,Vs);
            
            %Fix normal directions for caps if needed
            E1=patchEdges(Fs(Cs==1,:),0);
            for q=2:1:3
                Eb=patchBoundary(Fs(Cs==q,:),Vs);
                %If cap boundary edges are a member of the tube mesh edges the cap
                %faces are pointing the wrong way
                if any(isrowmember(Eb,E1,0))
                    Fs(Cs==q,:)=fliplr(Fs(Cs==q,:));
                end
            end
        end
end

%% Collect output

varargout{1}=Fs;
varargout{2}=Vs;
varargout{3}=Cs_C;
varargout{4}=Cs_rgb;
varargout{5}=Cs_d;
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
