function [varargout]=tet2hex(varargin)



%%
switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        tet2HexMethod=1; 
    case 3
        E=varargin{1};
        V=varargin{2};        
        tet2HexMethod=varargin{3};         
end

%% Mid edge sets
edgeMat=[E(:,[1 2]); E(:,[2 3]);  E(:,[3 1]);... %bottom         
         E(:,[1 4]); E(:,[2 4]);  E(:,[3 4]);];   %top-bottom edges

E_sort=sort(edgeMat,2); %Sorted edges matrix
[~,ind1,~]=unique(E_sort,'rows');
edgeMat=edgeMat(ind1,:);

numPoints = size(V,1);
numEdges = size(edgeMat,1);

% Get indices of the four edges associated with each face
A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
A = max(A,A'); %Copy symmetric

%Indices for A matrix
indA_12=E(:,1)+(E(:,2)-1)*numPoints;
indA_23=E(:,2)+(E(:,3)-1)*numPoints;
indA_31=E(:,3)+(E(:,1)-1)*numPoints;

indA_14=E(:,1)+(E(:,4)-1)*numPoints;
indA_24=E(:,2)+(E(:,4)-1)*numPoints;
indA_34=E(:,3)+(E(:,4)-1)*numPoints;

%Get indices for vertex array
indV_12=full(A(indA_12));
indV_23=full(A(indA_23));
indV_31=full(A(indA_31));

indV_14=full(A(indA_14));
indV_24=full(A(indA_24));
indV_34=full(A(indA_34));

%% Mid element
indV_mid=(1:1:size(E,1))'+numPoints+size(edgeMat,1);

%% Mid face

%Element faces matrix
F =[E(:,[1 2 3]);... %top
    E(:,[1 2 4]);... %side 1
    E(:,[2 3 4]);... %side 2
    E(:,[3 1 4]);... %side 3
    ]; 

F_sort=sort(F,2); %Sorted edges matrix
[~,ind1,ind2]=unique(F_sort,'rows');
F=F(ind1,:);

indV_midFace=(1:1:size(F_sort,1))';
indOffset=numPoints+size(edgeMat,1)+size(E,1);

indV_midFace123=ind2(indV_midFace((1-1)*size(E,1)+(1:size(E,1))))+indOffset;
indV_midFace124=ind2(indV_midFace((2-1)*size(E,1)+(1:size(E,1))))+indOffset;
indV_midFace234=ind2(indV_midFace((3-1)*size(E,1)+(1:size(E,1))))+indOffset;
indV_midFace314=ind2(indV_midFace((4-1)*size(E,1)+(1:size(E,1))))+indOffset;

%% Create element array

 Es=[E(:,1) indV_12 indV_midFace123 indV_31 indV_14  indV_midFace124 indV_mid indV_midFace314;...%Corner hex 1
    indV_12  E(:,2) indV_23 indV_midFace123 indV_midFace124 indV_24 indV_midFace234 indV_mid;...%Corner hex 2
     indV_midFace123 indV_23 E(:,3) indV_31 indV_mid indV_midFace234 indV_34 indV_midFace314;...%Corner hex 3
     indV_14  indV_midFace124 indV_mid indV_midFace314 E(:,4) indV_24 indV_midFace234 indV_34;...%Corner hex 4
    ];
% Es=Es(:,[4 3 2 1 8 7 6 5]);

%% Create vertex arrays

%new mid-edge points
Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:));
     
switch tet2HexMethod
    case 1        
        %new mid-element points
        Vm=zeros(size(E,1),3);
        for q=1:1:size(V,2)
            X=V(:,q);
            if size(E,1)==1
                Vm(:,q)=mean(X(E)',2);
            else
                Vm(:,q)=mean(X(E),2);
            end
        end
        
        %new mid-face points
        Vf=zeros(size(F,1),3);
        for q=1:1:size(V,2)
            X=V(:,q);
            Vf(:,q)=mean(X(F),2);
        end
    case 2        
        %new mid-element points
        TET=triangulation(E,V);
        Vm = incenter(TET,(1:size(E,1))');
        
        %new mid-face points
        TR=triangulation(F,V);
        Vf = incenter(TR,(1:size(F,1))');     
end

Vs=[V; Vn; Vm; Vf]; %Join point sets

CVs=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1); 3*ones(size(Vf,1),1);];

%%

varargout{1}=Es;
varargout{2}=Vs;
varargout{3}=CVs;

%% old
% %% Get faces
% F=[TET(:,[2 1 3]); TET(:,[1 2 4]); TET(:,[2 3 4]); TET(:,[3 1 4])];
% 
% %% Deriving new coordinate sets
% 
% %The original vertices
% X=V(:,1); Y=V(:,2); Z=V(:,3);
% 
% %The mid-face points
% Vmf=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];
% 
% %The mid points
% if size(TET,1)==1
%     Vm=[mean(X(TET)) mean(Y(TET)) mean(Z(TET))];
% else
%     Vm=[mean(X(TET),2) mean(Y(TET),2) mean(Z(TET),2)];
% end
% 
% %The mid edge points
% indEdge1=[1 2]; indEdge2=[2 3]; indEdge3=[3 1];
% indEdge4=[1 4]; indEdge5=[2 4]; indEdge6=[3 4];
% edgeIndAll=[TET(:,indEdge1);TET(:,indEdge2);TET(:,indEdge3);...
%     TET(:,indEdge4);TET(:,indEdge5);TET(:,indEdge6);];
% 
% if size(TET,1)==1
%     Xt=reshape(X(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
%     Yt=reshape(Y(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
%     Zt=reshape(Z(edgeIndAll),size(edgeIndAll,1),size(edgeIndAll,2));
% else
%     Xt=X(edgeIndAll);
%     Yt=Y(edgeIndAll);
%     Zt=Z(edgeIndAll);
% end
% Vme=[mean(Xt,2) mean(Yt,2) mean(Zt,2)];
% 
% %the joint coordinate set
% Vhex=[V; Vme; Vmf; Vm];
% 
% %% Creating hexahedral elements
% 
% %Numbers for fixing point indices
% numVert=size(V,1);
% numTET=size(TET,1);
% numVme=size(Vme,1);
% numVmf=size(Vmf,1);
% 
% %Indices for corner, edge-, face- points etc.
% mixInd=[1 1 1 3 4 2 0 4;...
%         2 2 1 1 5 3 0 2;...
%         3 3 1 2 6 4 0 3;...
%         4 4 2 5 6 4 0 3];
% 
% HEX=zeros(size(TET,1).*4,8); %The matrix for the hexahedral element spec.
% for hexInd=1:1:4
%     
%     %Corner
%     ind1=TET(:,mixInd(hexInd,1));
%     
%     %Mid edge point
%     edgeInd=mixInd(hexInd,2);
%     ind2=numVert+((edgeInd-1).*numTET)+(1:numTET);
%     
%     %Mid face point
%     faceInd=mixInd(hexInd,3);
%     ind3=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
%     
%     %Mid edge point
%     edgeInd=mixInd(hexInd,4);
%     ind4=numVert+((edgeInd-1).*numTET)+(1:numTET);
%     
%     %Mid edge point
%     edgeInd=mixInd(hexInd,5);
%     ind5=numVert+((edgeInd-1).*numTET)+(1:numTET);
%     
%     %Mid face point
%     faceInd=mixInd(hexInd,6);
%     ind6=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
%     
%     %Mid point
%     ind7=numVert+numVme+numVmf+(1:numTET);
%     
%     %Mid face point
%     faceInd=mixInd(hexInd,8);
%     ind8=numVert+numVme+((faceInd-1).*numTET)+(1:numTET);
%     
%     h=[ind1(:) ind2(:) ind3(:) ind4(:) ind5(:) ind6(:) ind7(:) ind8(:)];       
%     
%     if hexInd==4
%         h=h(:,[5 6 7 8 1 2 3 4]);
%     end
%     startInd=1+(hexInd-1).*(size(TET,1));
%     endInd=startInd-1+(size(TET,1));
%     HEX(startInd:endInd,:)=h;
%     
% end

% %% Removing double VERTICES
% numKeep=6;
% 
% [~,ind1,ind2]=unique(pround(Vhex,numKeep),'rows');
% 
% Vhex=Vhex(ind1,:);
% HEX=ind2(HEX);
 
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
