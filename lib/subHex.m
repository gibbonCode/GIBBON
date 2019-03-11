function [varargout]=subHex(varargin)

% function [E,V,C,CV]=subHex(E,V,n,splitMethod)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
%
% 2017/08/31 Added function description to top of function
% ------------------------------------------------------------------------


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
if n>0
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
                Vm=patchCentre(E,V); %Mid element nodes
                
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
                Vs = [V; Vn;]; %Join point sets
                CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1);];
                
        end
        
        %% Override input
        C=repmat(C,[size(Es,1)/size(E,1),1]);
        E=Es;
        V=Vs;
    end
end

varargout{1}=E;
varargout{2}=V;
varargout{3}=C;
varargout{4}=CV;


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
