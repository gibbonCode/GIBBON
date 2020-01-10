function [outputStruct]=subTriLocal(inputStruct)

% function [outputStruct]=subTriLocal(inputStruct)
% ------------------------------------------------------------------------
% This function refines (using the subTri function) the triangles in F
% defined by indFaces and splits neighbouring triangles to create a
% conforming mesh. This way local refinement can be achieved. 
% The output mesh is stored in Fq and Vq and a face color list Cq can also
% be requested which lists whether a triangle is unaltered (Cq==1), is
% subdevided (Cq==2) or has been split to connect the two regions (Cq==3).
% The optional input f (default is 0) defines the location of the new
% points introduced for the transition elements. Using f>0 (and <1) will
% place these points closer to the coarse mesh nodes. The optional output
% indInitial is a list containing all the original nodes. 
%
%
% See also: subTri
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/12/17 Created
%------------------------------------------------------------------------

%% Parse input

if ~isfield(inputStruct,'F')
    error('Missing field in input structure: F');
else
    F=inputStruct.F; 
end

if ~isfield(inputStruct,'V')
    error('Missing field in input structure: V');
else
    V=inputStruct.V; 
end

if ~isfield(inputStruct,'indFaces')
    error('Missing field in input structure: indFaces');
else
    indFaces=inputStruct.indFaces; 
end

if ~isfield(inputStruct,'f')
    f=0;
else
    f=inputStruct.f; 
end

if ~isfield(inputStruct,'C')
    C=[];
else
    C=inputStruct.C; 
end

% Check validity of f
if f>=1 || f<0 
   error('The factor f should f>=0 and f<1'); 
end

%% Check for enclosed triangles and add them to the list

logicFaces=false(size(F,1),1);
logicFaces(indFaces)=1;
Fs=F(logicFaces,:);

ind_Fs=unique(Fs(:));
logicFaces=sum(ismember(F,ind_Fs),2)==3;

%% The refined/subdevided face set

Fs=F(logicFaces,:);
ind_Fs=unique(Fs(:));

[Fs,Vs,~]=patchCleanUnused(Fs,V);
numPoints=size(Vs,1);
[Fs,Vs]=subtri(Fs,Vs,1,1);

%Treat color bookkeeping
if ~isempty(C)
    Cs=C(logicFaces,:);
    Cs=repmat(Cs,[4,1]);
end

%% The remaining faces

Fr=F(~logicFaces,:);
logicEdgeFaces=sum(ismember(Fr,ind_Fs),2)==2; 

%Treat color bookkeeping
if ~isempty(C)
    Cr=C(~logicFaces,:);
end

%% The faces which will be kept
Fk=Fr(~logicEdgeFaces,:);
[Fk,Vk,~]=patchCleanUnused(Fk,V);

%Treat color bookkeeping
if ~isempty(C)
    Ck=Cr(~logicEdgeFaces,:);
end

%% The faces that will be cut
Fc=Fr(logicEdgeFaces,:);
ind_Fc=unique(Fc(:));
indEdgePoints=ind_Fc(ismember(ind_Fc,ind_Fs)); 

[Fc,Vc,indFix2]=patchCleanUnused(Fc,V);
indEdgePoints=indFix2(indEdgePoints);
logicEdge=ismember(Fc,indEdgePoints);

%Treat color bookkeeping
if ~isempty(C)
    Cc=Cr(logicEdgeFaces,:);
end

%Compose new point set
X=Vc(:,1);
XF=X(Fc).*logicEdge;
XF(~logicEdge)=NaN;
Xn=gnanmean(XF,2);
Y=Vc(:,2);
YF=Y(Fc).*logicEdge;
YF(~logicEdge)=NaN;
Yn=gnanmean(YF,2);
Z=Vc(:,3);
ZF=Z(Fc).*logicEdge;
ZF(~logicEdge)=NaN;
Zn=gnanmean(ZF,2);
Vn=[Xn Yn Zn];

Vct=[Vc;Vn]; %The total node set for cut faces

ind_Vn=(1:size(Vn,1))+size(Vc,1);

logic12=logicEdge(:,1) & logicEdge(:,2); 
F12=Fc(logic12,:);
ind12=find(logic12)+size(Vc,1);
F12n=[F12(:,1) ind12 F12(:,3);ind12 F12(:,2)  F12(:,3)];

logic13=logicEdge(:,1) & logicEdge(:,3); 
F13=Fc(logic13,:);
ind13=find(logic13)+size(Vc,1);
F13n=[F13(:,1) F13(:,2) ind13; F13(:,2) F13(:,3) ind13;];

logic23=logicEdge(:,2) & logicEdge(:,3); 
F23=Fc(logic23,:);
ind23=find(logic23)+size(Vc,1);
F23n=[F23(:,1) F23(:,2) ind23; F23(:,1) ind23 F23(:,3);];

Fct=[F12n;F13n;F23n]; %The faces for the cut face set

%Treat color bookkeeping
if ~isempty(C)
    Cct=[Cc(logic12,:); Cc(logic12,:); ...
         Cc(logic13,:); Cc(logic13,:); ...
         Cc(logic23,:); Cc(logic23,:);];
end

%% Merge meshes

Vq=[Vk;Vs;Vct];
Fq=[Fk;Fs+size(Vk,1);Fct+size(Vk,1)+size(Vs,1)];
if ~isempty(C)
    Cq=[Ck;Cs;Cct];
end
faceTypeLabel=[ones(size(Fk,1),1); 2*ones(size(Fs,1),1); 3*ones(size(Fct,1),1)];
ind_Vn=ind_Vn+size(Vk,1)+size(Vs,1);
indInitial=[1:size(Vk,1) (1:numPoints)+size(Vk,1)];

%Merge non-unique nodes
[Fq,Vq,~,ind2]=mergeVertices(Fq,Vq);
ind_Vn=ind2(ind_Vn);
indInitial=ind2(indInitial);

%% Displace point using f

if f>0
    Fw3=Fq(faceTypeLabel==3,:);
    E3=patchEdges(Fw3,1);
    L3=ismember(E3,ind_Vn);
    
    Fw2=Fq(faceTypeLabel==2,:);
    E2=patchEdges(Fw2,1);
    L2=ismember(E3,E2(:));
    
    LE=any(L3,2) & any(~L2,2);
    LE3=L3(LE,:);
    E=E3(LE,:);
    [~,indSort]=sort(LE3,2);
    
    ind1=find(indSort==1);
    ind2=find(indSort==2);
    [ind1Rows,~] = ind2sub(size(E),ind1);
    [ind2Rows,~] = ind2sub(size(E),ind2);
    
    En=E;
    En(ind1Rows,1)=E(ind1);
    En(ind2Rows,2)=E(ind2);
    
    X=Vq(:,1);
    Y=Vq(:,2);
    Z=Vq(:,3);
    
    XE=X(En);
    YE=Y(En);
    ZE=Z(En);

    XE=(1-f).*XE(:,2)+f*XE(:,1);
    YE=(1-f).*YE(:,2)+f*YE(:,1);
    ZE=(1-f).*ZE(:,2)+f*ZE(:,1);
    
    VE=[XE YE ZE];
    Vq(En(:,2),:)=VE;
end

%% Smooth
% 
% if ~isempty(smoothPar)
%     smoothPar.RigidConstraints=indInitial;    
%     [Vq]=tesSmooth(Fq,Vq,[],smoothPar);
% end

%%

outputStruct.F=Fq;
outputStruct.V=Vq;
if ~isempty(C)
    outputStruct.C=Cq;
end
outputStruct.faceTypeLabel=faceTypeLabel;
outputStruct.indInitial=indInitial;
 
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
