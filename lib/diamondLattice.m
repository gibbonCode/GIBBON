function [Ep,Et,VT,Ct]=diamondLattice(sampleSize,nRepeat,strutThickness,plotOn)


%%



%% Derive shrinkfactor

shrinkFactor=strutThickness./((sampleSize./nRepeat).*(sqrt(2)./2));

if shrinkFactor>1
    error('shrinkFactor cannot exceed 1');
end

%% Construct an octahedron

%Vertices of octahedron
V=[-1  0  0;...
    0 -1  0;...
    1  0  0;...
    0  1  0;...
    0  0 -1;...
    0  0  1]/2;

%Faces of octahedron
F=[5 2 1; 5 3 2; 5 4 3; 5 1 4; 1 2 6; 2 3 6; 3 4 6; 4 1 6;];

%% Create tetrahedra on faces of octahedron
NS=patchNormal(F,V).*sqrt(6)/3./sqrt(2);
VT=patchCentre(F,V);
VS=NS+VT;
Vs=[V;VS];
Es=[F (1:size(F,1))'+size(V,1)];
[Fs]=element2patch(Es);

%%
if plotOn==1
    cFigure; hold on;
    gpatch(F,V,'bw','b',1);
    gpatch(Fs,Vs,'rw','r',0.25);
    axisGeom;
    camlight headlight;
    drawnow;
end

%% Offset octahedra

nCopies=nRepeat*ones(1,3)+1; %Number of "copies" in each direction

V_offsets=[1 0 0; 0 1 0; 0 0 1]; %Offset vectors
nCopies1=nCopies;
[Vg1]=gridClone(V,nCopies1,V_offsets);

nCopies2=nCopies;
nCopies2(1:2)=nCopies2(1:2)-1;
[Vg2]=gridClone(V,nCopies2,V_offsets);
 
nCopies3=nCopies;
nCopies3(1)=nCopies3(1)-1;
nCopies3(3)=nCopies3(3)-1;
[Vg3]=gridClone(V,nCopies3,V_offsets);

nCopies4=nCopies;
nCopies4(2)=nCopies4(2)-1;
nCopies4(3)=nCopies4(3)-1;
[Vg4]=gridClone(V,nCopies4,V_offsets);

Vg2(:,1)=Vg2(:,1)+0.5;
Vg2(:,2)=Vg2(:,2)+0.5;

Vg3(:,1)=Vg3(:,1)+0.5;
Vg3(:,3)=Vg3(:,3)+0.5;

Vg4(:,2)=Vg4(:,2)+0.5;
Vg4(:,3)=Vg4(:,3)+0.5;

V1=[Vg1;Vg2;Vg3;Vg4;];
 
E1=reshape((1:1:size(V1,1)),6,size(V1,1)/6)';
C1=(1:1:size(E1,1))'; %Element colors
[F1,CF1]=element2patch(E1,C1,'octa6');
V_centre1=patchCentre(E1,V1);

%% Offset tetrahedra

V_offsets=[1 0 0; 0 1 0; 0 0 1]; %Offset vectors

V2=gridClone(Vs,nCopies,V_offsets);

n=prod(nCopies);
E2=repmat(Es,n,1);

ind=((0:1:n-1)'*ones(1,8))';

ind=ind(:).*size(Vs,1);
E2=E2+ind(:,ones(1,4));

C2=(1:1:size(E2,1))'; %Element colors

F2=[E2(:,[1 2 3]);... %face 1 2 3
    E2(:,[1 2 4]);... %face 1 2 4
    E2(:,[2 3 4]);... %face 2 3 4
    E2(:,[3 1 4])];   %face 1 3 4

CF2=repmat(C2,4,1); %Replicate color data

V_centre2=patchCentre(E2,V2);

%%

if plotOn==1
cFigure; 
subplot(1,2,1); hold on; 
title('Repeated octahedra')
gpatch(F1,V1,CF1,'k',0.8);
plotV(V_centre1,'k.','MarkerSize',25)
colormap gjet; 
axisGeom; 
camlight headlight; 

subplot(1,2,2); hold on; 
title('Repeated tetrahedra')
gpatch(F2,V2,CF2,'k',0.8);
plotV(V_centre2,'k.','MarkerSize',25)
colormap gjet; 
axisGeom; 
camlight headlight; 

drawnow;
end

%%

[E1,V1,V1s]=scalePatch(E1,V1,shrinkFactor);
[E2,V2,V2s]=scalePatch(E2,V2,shrinkFactor);

F1=element2patch(E1,[],'octa6');
F2=element2patch(E2,[],'tet4');

VTs=[V1s;V2s]; %Join point sets
VT=[V1;V2]; %Join point sets
F12=[F1;F2+size(V1,1)];
E2=E2+size(V1,1); %Offset indices of tetrahedral set
F2=F2+size(V1,1); %Offset indices of tetrahedral set
[F12s,~]=mergeVertices(F12,VTs); %Merge vertices

F12_sort=sort(F12s,2);
[~,indUni,ind2]=unique(F12_sort,'rows');
I1=ind2(1:size(F1,1));
I2=ind2(size(F1,1)+1:end);
ind_I1=1:1:numel(I1);
ind_I2=1:1:numel(I2);
A=zeros(size(indUni));
A(I1)=ind_I1;
B=zeros(size(indUni));
B(I2)=ind_I2;
C=[A(:) B(:)];
logicPair=all(C>0,2);
C=C(logicPair,:);
CC=[I1;I2];

%%
if plotOn==1
cFigure; 
gpatch(F1,VT,'bw','b',0.25);
gpatch(F2,VT,'rw','r',0.25);
colormap gjet; 
axisGeom; 
camlight headlight; 
drawnow;
end
%%

F1p=F1(C(:,1),:);
F2p=F2(C(:,2),:);

forder=[1 2 3; 2 3 1; 3 1 2; 3 2 1; 2 1 3; 1 3 2;];
D=zeros(size(F1p,1),size(forder,1));
for q1=1:1:size(forder,1)
    f=forder(q1,:);
    d=zeros(size(F1p,1),1);
    for q=1:1:3
        X=VT(:,q);
        d=d+(X(F1p)-X(F2p(:,f))).^2;
    end
    d=sum(sqrt(d),2);
    D(:,q1)=d;    
    
end

[~,indMin]=min(D,[],2);

for q=1:1:max(indMin)
    logicNow=indMin==q;
    F2p(logicNow,:)=F2p(logicNow,forder(q,:));
end

Ep=[F1p F2p];

%%

[E1t,Vt]=octa2tet(E1,VT);
E1t=E1t+size(VT,1);
VT=[VT;Vt];
[F12,VT,~,ind2]=mergeVertices(F12,VT);
E1t=ind2(E1t);
E2=ind2(E2);
Ep=ind2(Ep);

%%
if plotOn==1
    Fp=element2patch(Ep,[],'penta6');
    F1t=element2patch(E1t,[],'tet4');
    F2=element2patch(E2,[],'tet4');
    
    cFigure; hold on;
    gpatch(Fp,VT,'rw','r',0.5);
    gpatch(F1t,VT,'gw','g',0.5);
    gpatch(F2,VT,'bw','b',0.5);
    
    axisGeom;
    camlight headlight;
    drawnow;
end

%%
VE1=patchCentre(E1t,VT);
logicKeep_E1t=all(VE1>0 & VE1<nRepeat,2);
E1t=E1t(logicKeep_E1t,:);

VE2=patchCentre(E2,VT);
logicKeep_E2=all(VE2>0 & VE2<nRepeat,2);
E2=E2(logicKeep_E2,:);

VEp=patchCentre(Ep,VT);
logicKeep_Ep=all(VEp>0 & VEp<nRepeat,2);
Ep=Ep(logicKeep_Ep,:);

indUsed=unique([Ep(:);E1t(:);E2(:)]);
[~,VT,indFix]=patchCleanUnused(indUsed,VT); 
E1t=indFix(E1t);
E2=indFix(E2);
Ep=indFix(Ep);

Fp=element2patch(Ep,[],'penta6');
F1t=element2patch(E1t,[],'tet4');
F2=element2patch(E2,[],'tet4');

if plotOn==1
    cFigure; hold on;
    gpatch(Fp,VT,'rw','r',0.5);
    gpatch(F1t,VT,'gw','g',0.5);
    gpatch(F2,VT,'bw','b',0.5);
    
    axisGeom;
    camlight headlight;
    drawnow;
end

%% Scaling model to desired size

VT=VT./nRepeat;
VT=VT.*sampleSize;

%%

Et=[E1t;E2];
Ct=[ones(size(E1t,1),1); 2*ones(size(E2,1),1);];

end


function [Vg]=gridClone(Vs,nCopies,V_offsets)

nTotal=prod(nCopies); %Total number of copies
indC=ones(size(Vs,1),1)*(1:1:nTotal);
indC=indC(:);
[I,J,K] = ind2sub(nCopies,indC);
I=I-1; J=J-1; K=K-1;

%Defining offsets
D1=I*V_offsets(1,:);
D2=J*V_offsets(2,:);
D3=K*V_offsets(3,:);

%Defining vertices matrix
Vg=repmat(Vs,nTotal,1)+(D1+D2+D3);

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
