function [varargout]=subHex(varargin)


%% Parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        n=1;
        splitMethod=1;
    case 3
        E=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=1;
    case 4
        E=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=varargin{4};
end

C=(1:1:size(E,1))'; %Element colors or indices
CV=[];

%%
for qIter=1:1:n
    switch splitMethod
        case 1
            %% Mid edge sets
            edgeMat=[E(:,[1 2]); E(:,[2 3]);  E(:,[3 4]); E(:,[4 1]);... %top
                E(:,[5 6]); E(:,[6 7]);  E(:,[7 8]); E(:,[8 5]);... %bottom
                E(:,[1 5]); E(:,[2 6]);  E(:,[3 7]); E(:,[4 8])];   %top-bottom edges
            
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
            indA_34=E(:,3)+(E(:,4)-1)*numPoints;
            indA_41=E(:,4)+(E(:,1)-1)*numPoints;
            
            indA_56=E(:,5)+(E(:,6)-1)*numPoints;
            indA_67=E(:,6)+(E(:,7)-1)*numPoints;
            indA_78=E(:,7)+(E(:,8)-1)*numPoints;
            indA_85=E(:,8)+(E(:,5)-1)*numPoints;
            
            indA_15=E(:,1)+(E(:,5)-1)*numPoints;
            indA_26=E(:,2)+(E(:,6)-1)*numPoints;
            indA_37=E(:,3)+(E(:,7)-1)*numPoints;
            indA_48=E(:,4)+(E(:,8)-1)*numPoints;
            
            %Get indices for vertex array
            indV_12=full(A(indA_12));
            indV_23=full(A(indA_23));
            indV_34=full(A(indA_34));
            indV_41=full(A(indA_41));
            
            indV_56=full(A(indA_56));
            indV_67=full(A(indA_67));
            indV_78=full(A(indA_78));
            indV_85=full(A(indA_85));
            
            indV_15=full(A(indA_15));
            indV_26=full(A(indA_26));
            indV_37=full(A(indA_37));
            indV_48=full(A(indA_48));
            
            %% Mid element
            indV_mid=(1:1:size(E,1))'+numPoints+size(edgeMat,1);
            
            %% Mid face
            
            %Element faces matrix
            F =[E(:,[4 3 2 1]);... %top
                E(:,[5 6 7 8]);... %bottom
                E(:,[1 2 6 5]);... %side 1
                E(:,[3 4 8 7]);... %side 2
                E(:,[2 3 7 6]);... %front
                E(:,[5 8 4 1]);]; %back
            
            F_sort=sort(F,2); %Sorted edges matrix
            [~,ind1,ind2]=unique(F_sort,'rows');
            F=F(ind1,:);
            
            indV_midFace=(1:1:size(F_sort,1))';
            indOffset=numPoints+size(edgeMat,1)+size(E,1);
            
            indV_midFace4321=ind2(indV_midFace((1-1)*size(E,1)+(1:size(E,1))))+indOffset;
            indV_midFace5678=ind2(indV_midFace((2-1)*size(E,1)+(1:size(E,1))))+indOffset;
            indV_midFace1265=ind2(indV_midFace((3-1)*size(E,1)+(1:size(E,1))))+indOffset;
            indV_midFace3487=ind2(indV_midFace((4-1)*size(E,1)+(1:size(E,1))))+indOffset;
            indV_midFace2376=ind2(indV_midFace((5-1)*size(E,1)+(1:size(E,1))))+indOffset;
            indV_midFace5841=ind2(indV_midFace((6-1)*size(E,1)+(1:size(E,1))))+indOffset;
            
            %% Create element array
            Es=[E(:,1) indV_12 indV_midFace4321 indV_41 indV_15 indV_midFace1265 indV_mid indV_midFace5841;...%Corner hex 1
                indV_12 E(:,2) indV_23 indV_midFace4321 indV_midFace1265 indV_26 indV_midFace2376 indV_mid;... %Corner hex 2
                indV_midFace4321 indV_23 E(:,3) indV_34 indV_mid indV_midFace2376 indV_37 indV_midFace3487;... %Corner hex 3
                indV_41 indV_midFace4321 indV_34 E(:,4) indV_midFace5841 indV_mid indV_midFace3487 indV_48;... %Corner hex 4
                indV_15 indV_midFace1265 indV_mid indV_midFace5841 E(:,5) indV_56 indV_midFace5678 indV_85 ;...%Corner hex 5
                indV_midFace1265 indV_26 indV_midFace2376 indV_mid  indV_56 E(:,6) indV_67 indV_midFace5678;... %Corner hex 6
                indV_mid indV_midFace2376 indV_37 indV_midFace3487 indV_midFace5678 indV_67 E(:,7) indV_78 ;... %Corner hex 7
                indV_midFace5841 indV_mid indV_midFace3487 indV_48 indV_85 indV_midFace5678 indV_78 E(:,8) ;... %Corner hex 8
                ];
            
            %% Create vertex array
            Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
            Vm=zeros(size(E,1),3);
            for q=1:1:size(V,2)
                X=V(:,q);
                if size(E,1)==1
                    Vm(:,q)=mean(X(E)',2);
                else
                    Vm(:,q)=mean(X(E),2);
                end
            end
            Vf=zeros(size(F,1),3);
            for q=1:1:size(V,2)
                X=V(:,q);
                Vf(:,q)=mean(X(F),2);
            end
            Vs = [V; Vn; Vm; Vf]; %Join point sets
            
            CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1); 3*ones(size(Vf,1),1);];

        case 2
            
            %% Mid edge sets
            edgeMat=[E(:,[1 5]); E(:,[2 6]);  E(:,[3 7]); E(:,[4 8])];   %top-bottom edges
            
            E_sort=sort(edgeMat,2); %Sorted edges matrix
            [~,ind1,~]=unique(E_sort,'rows');
            edgeMat=edgeMat(ind1,:);
            
            numPoints = size(V,1);
            numEdges = size(edgeMat,1);
            
            % Get indices of the four edges associated with each face
            A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
            A = max(A,A'); %Copy symmetric
            
            %Indices for A matrix
            indA_15=E(:,1)+(E(:,5)-1)*numPoints;
            indA_26=E(:,2)+(E(:,6)-1)*numPoints;
            indA_37=E(:,3)+(E(:,7)-1)*numPoints;
            indA_48=E(:,4)+(E(:,8)-1)*numPoints;
            
            %Get indices for vertex array
            indV_15=full(A(indA_15));
            indV_26=full(A(indA_26));
            indV_37=full(A(indA_37));
            indV_48=full(A(indA_48));
            
            %% Create element array
            Es=[E(:,1:4) indV_15 indV_26 indV_37 indV_48;...
                indV_15 indV_26 indV_37 indV_48 E(:,5:8)];
            
            %% Create vertex array
            Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
            Vm=zeros(size(E,1),3);
            for q=1:1:size(V,2)
                X=V(:,q);
                if size(E,1)==1
                    Vm(:,q)=mean(X(E)',2);
                else
                    Vm(:,q)=mean(X(E),2);
                end
            end
            
            Vs = [V; Vn;]; %Join point sets
            
            CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1);];
            
        case 3
            
            %% Mid edge sets
            edgeMat=[E(:,[1 2]); E(:,[5 6]);  E(:,[7 8]); E(:,[3 4])]; %left-right edges
            
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
            indA_56=E(:,5)+(E(:,6)-1)*numPoints;
            indA_78=E(:,7)+(E(:,8)-1)*numPoints;
            indA_34=E(:,3)+(E(:,4)-1)*numPoints;
            
            %Get indices for vertex array
            indV_12=full(A(indA_12));
            indV_56=full(A(indA_56));
            indV_78=full(A(indA_78));
            indV_34=full(A(indA_34));
            
            %% Create element array
            Es=[E(:,1) indV_12 indV_34 E(:,4) E(:,5) indV_56 indV_78 E(:,8);...
                indV_12 E(:,2) E(:,3) indV_34 indV_56  E(:,6) E(:,7) indV_78];
            
            %% Create vertex array
            Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
            Vm=zeros(size(E,1),3);
            for q=1:1:size(V,2)
                X=V(:,q);
                if size(E,1)==1
                    Vm(:,q)=mean(X(E)',2);
                else
                    Vm(:,q)=mean(X(E),2);
                end
            end
            
            Vs = [V; Vn;]; %Join point sets
            CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1);];
            
        case 4
            
            %% Mid edge sets
            edgeMat=[E(:,[4 1]); E(:,[8 5]);  E(:,[6 7]); E(:,[2 3])]; %front-back edges
            
            E_sort=sort(edgeMat,2); %Sorted edges matrix
            [~,ind1,~]=unique(E_sort,'rows');
            edgeMat=edgeMat(ind1,:);
            
            numPoints = size(V,1);
            numEdges = size(edgeMat,1);
            
            % Get indices of the four edges associated with each face
            A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
            A = max(A,A'); %Copy symmetric
            
            %Indices for A matrix
            indA_41=E(:,4)+(E(:,1)-1)*numPoints;
            indA_85=E(:,8)+(E(:,5)-1)*numPoints;
            indA_67=E(:,6)+(E(:,7)-1)*numPoints;
            indA_23=E(:,2)+(E(:,3)-1)*numPoints;
            
            %Get indices for vertex array
            indV_41=full(A(indA_41));
            indV_85=full(A(indA_85));
            indV_67=full(A(indA_67));
            indV_23=full(A(indA_23));                       
            
            %% Create element array
            Es=[E(:,[1 2]) indV_23 indV_41 E(:,[5 6]) indV_67 indV_85;...
                indV_41 indV_23 E(:,[3 4]) indV_85 indV_67 E(:,[7 8])];
            
            %% Create vertex array
            Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
            Vm=zeros(size(E,1),3);
            for q=1:1:size(V,2)
                X=V(:,q);
                if size(E,1)==1
                    Vm(:,q)=mean(X(E)',2);
                else
                    Vm(:,q)=mean(X(E),2);
                end
            end
            
            Vs = [V; Vn;]; %Join point sets
            CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1);];
            
    end
        
    C=repmat(C,[size(Es,1)/size(E,1),1]);
    E=Es;
    V=Vs;    
end

varargout{1}=E; 
varargout{2}=V; 
varargout{3}=C; 
varargout{4}=CV; 

 
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
