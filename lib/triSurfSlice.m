function [varargout]=triSurfSlice(varargin)
% function [F,V,C,logicSide,Eb]=triSurfSlice(F,V,C,P,n,snapTolerance)
%-------------------------------------------------------------------------
% This function slices the input mesh (defined by faces F, vertices V and
% optional color data C) using a plane specified by a point P and a normal
% vector to the plane n. Edges crossing the plane are cut at the plane, and
% new points are introduced at the intersection, and new triangles are
% created. Points within snapTolerance from the plane are snapped to the
% plane.
%
% 2020/04/10 Created
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=ones(size(F,1),1);
        P=mean(V,1);
        n=[0 0 1];
        snapTolerance=[];
        logicExclude=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        P=mean(V,1);
        n=[0 0 1];
        snapTolerance=[];
        logicExclude=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        P=varargin{4};
        n=[0 0 1];
        snapTolerance=[];
        logicExclude=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        P=varargin{4};
        n=varargin{5};
        snapTolerance=[];
        logicExclude=[];
    case 6
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        P=varargin{4};
        n=varargin{5};
        snapTolerance=varargin{6};
        logicExclude=[];
    case 7
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
        P=varargin{4};
        n=varargin{5};
        snapTolerance=varargin{6};
        logicExclude=varargin{7};
end

%Check color
if isempty(C)
    C=ones(size(F,1),1);
end

%Check P
if isempty(P)
    P=mean(V,1); %Use default as mean of coordinates
end

%Check normal direction
if isempty(n)
    n=[0 0 1]; %Default z-direction
end
n=vecnormalize(n); %Normalize n vector

%Use default for snap tolerance if empty
if isempty(snapTolerance)
    %Default 1/100 of mean edge length
    snapTolerance=mean(patchEdgeLengths(F,V))/100; 
end

if ~isempty(logicExclude)
    indExclude=unique(F(logicExclude,:));
else
    indExclude=[];
end

%% Construct rotation matrix
nz=[0 0 1]; %The z-vector
R=vecPair2Rot(n,nz);

%% Derive connectivity arrays
[CC]=patchConnectivity(F,V,{'fe','ev'});
E=CC.edge.vertex; %The edges-vertex connectivity array
FE=CC.face.edge; %The face-edge connectivity array 

%%

%Rotate coordinates
Vr=V-P(ones(size(V,1),1),:);
Vr=Vr*R';

%Snap close points to plane
z=Vr(:,3);
logicSnap=abs(z)<=snapTolerance;

%Simple snap (distorts shape for large snap tolerance)
Vr(logicSnap,3)=0; 

% WIP: Alternative using push allong edges. Suffers from 0 dot product
% issues. 
% [~,~,NV]=patchNormal(F,V);
% nv=vecnormalize(NV(logicSnap,:)); %Mesh normal vectors
% nt=vecnormalize(cross(nv,nz(ones(nnz(logicSnap),1),:),2)); %Mesh tangent vectors
% nz_mesh=vecnormalize(cross(nt,nv,2)); %Mesh tangent close to z vec
% d=dot(nz_mesh,nz(ones(nnz(logicSnap),1),:),2);
% logicInvalid=abs(d)<eps(1);
% pushFactor=-Vr(logicSnap,3)./d;
% pushFactor(logicInvalid)=0;
% Vr(logicSnap,:)=Vr(logicSnap,:)+pushFactor.*nz_mesh; %Push to plane
% indSnap=find(logicSnap);
% Vr(indSnap(logicInvalid),3)=0;

%Create logic for points above the plane
logicAbove=Vr(:,3)>0;
logicAbove(logicSnap)=1;

%Sample logic for edge points
logicAboveEdge=logicAbove(E);
logicSnapEdge=any(logicSnap(E),2);

%Find edges that cross plane
logicCrossing=sum(logicAboveEdge,2)==1;
logicCrossing(logicSnapEdge)=0;

if ~isempty(indExclude)
    logicCrossing=logicCrossing & ~any(ismember(E,indExclude),2);    
end

indCrossing=find(logicCrossing);
E_cross=E(logicCrossing,:); %Keep only relevant edges
logicAboveEdge=logicAboveEdge(logicCrossing,:);

%Sort edges array so the points above are last
[~,J_sort]=sort(logicAboveEdge,2);
i=(1:1:size(E_cross,1))';
indSort=sub2ind(size(E_cross),[i i],J_sort);
E_cross=E_cross(indSort);

NE=vecnormalize(Vr(E_cross(:,2),:)-Vr(E_cross(:,1),:));

d=dot(NE,nz(ones(size(E_cross,1),1),:),2);
pushFactor=-Vr(E_cross(:,1),3)./d;
logicInvalid=d<eps(0);
pushFactor(logicInvalid)=0;
VN=Vr(E_cross(:,1),:)+pushFactor.*NE; %Push to plane
VN(logicInvalid,3)=0;

FE(~ismember(FE,indCrossing))=0;
indFix=zeros(size(E,1),1);
indFix(indCrossing)=1:1:numel(indCrossing);
FE(FE>0)=indFix(FE(FE>0));
logicValid=FE>0;
logic1=sum(double(logicValid),2)==1;
logicKeep=~any(logicValid,2);

logic11=logic1 & logicValid(:,1);
logic12=logic1 & logicValid(:,2);
logic13=logic1 & logicValid(:,3);

logic21=~logic1 & logicValid(:,1) & logicValid(:,3);
logic22=~logic1 & logicValid(:,1) & logicValid(:,2);
logic23=~logic1 & logicValid(:,2) & logicValid(:,3);

FE_V=FE;
FE_V(FE_V>0)=FE_V(FE_V>0)+size(V,1);
V=[Vr;VN];
Fn_tri=[F(logic11,1)    FE_V(logic11,1) F(logic11,3)   ;...
    FE_V(logic11,1) F(logic11,2)    F(logic11,3)   ;...
    F(logic12,1)    F(logic12,2)    FE_V(logic12,2);...
    FE_V(logic12,2) F(logic12,3)    F(logic12,1)   ;...
    F(logic13,1)    F(logic13,2)    FE_V(logic13,3);...
    F(logic13,2)    F(logic13,3)    FE_V(logic13,3);...
    F(logic21,1)    FE_V(logic21,1) FE_V(logic21,3);...
    F(logic22,2)    FE_V(logic22,2) FE_V(logic22,1);...
    F(logic23,3)    FE_V(logic23,3) FE_V(logic23,2);...
    ];
Cn_tri=[repmat(C(logic11,:),2,1);repmat(C(logic12,:),2,1); repmat(C(logic13,:),2,1);...
    C(logic21,:); C(logic22,:); C(logic23,:);];

Fn_quad=[FE_V(logic21,1) F(logic21,2)    F(logic21,3)    FE_V(logic21,3);...
    F(logic22,1)    FE_V(logic22,1) FE_V(logic22,2) F(logic22,3);...
    F(logic23,1)    F(logic23,2)    FE_V(logic23,2) FE_V(logic23,3);...
    ];
Cn_quad=[C(logic21,:);C(logic22,:);C(logic23,:)];

[Fn_tri_quad,~,Cn_tri_quad]=quad2tri(Fn_quad,V,'a',Cn_quad);

Fn=[Fn_tri;Fn_tri_quad];
Cn=[Cn_tri;Cn_tri_quad];
F=[F(logicKeep,:);Fn];
C=[C(logicKeep,:);Cn];

VF=patchCentre(F,V);
logicSide=VF(:,3)<0;

if ~isempty(logicExclude)
    logicExclude=[logicExclude(logicKeep,:); false(size(Fn,1),1)];
    logicSide(logicExclude)=0;
end

V=V*R; %Rotate back
V=V+P(ones(size(V,1),1),:); %Shift back

EB=patchBoundary(F); %The full set of boundary edges
Eb_cut=patchBoundary(F(logicSide,:)); %Set containing old+new boundary edges
Eb=Eb_cut(~all(ismember(Eb_cut,EB),2),:);

%% Gather output
varargout{1}=F; 
varargout{2}=V; 
varargout{3}=C; 
varargout{4}=logicSide; 
varargout{5}=Eb; 

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
