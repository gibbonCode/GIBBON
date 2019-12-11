function [varargout]=subQuad(varargin)

% function [Fs,Vs,C,CV]=subQuad(F,V,n,splitMethod)
% ------------------------------------------------------------------------
% Sub-devides the quadrilateral faces defined by the patch format data F
% (faces) and V (vertices). Each face is split n times using the specified
% split method (splitMethod).
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2010/06/01 Created
% 2017/11/29 Fixed single face input related bug 
% 2018/11/06 Added splitMethod for splitting in a certain direction
% 2018/11/06 Added varargin based input handling
% 2018/11/06 Added color data handling
% 2018/11/06 Added variable output handling
% ------------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        n=1;
        splitMethod=1;
    case 3
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=1;
    case 4
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        splitMethod=varargin{4};
end

C=(1:1:size(F,1))'; %Face colors or indices
CV=[];

%%

%Treat 2D case
if size(V,2)==2
    V(:,3)=0;
    crop2Dopt=1;
else
    crop2Dopt=0;
end

if n>0
    for qIter=1:1:n
        switch splitMethod
            case 1
                
                edgeMat=[F(:,[1 2]); F(:,[2 3]);  F(:,[3 4]); F(:,[4 1])]; %Edges matrix
                Es=sort(edgeMat,2); %Sorted edges matrix
                [~,ind1,~]=unique(Es,'rows');
                edgeMat=edgeMat(ind1,:);
                
                numPoints = size(V,1);
                numEdges = size(edgeMat,1);
                
                % Get indices of the three edges associated with each face
                A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
                A = max(A,A'); %Copy symmetric
                
                %Indices for A matrix
                indA_12=F(:,1)+(F(:,2)-1)*numPoints;
                indA_23=F(:,2)+(F(:,3)-1)*numPoints;
                indA_34=F(:,3)+(F(:,4)-1)*numPoints;
                indA_41=F(:,4)+(F(:,1)-1)*numPoints;
                
                %Get indices for vertex array
                indV_12=full(A(indA_12));
                indV_23=full(A(indA_23));
                indV_34=full(A(indA_34));
                indV_41=full(A(indA_41));
                
                indV_mid=(1:1:size(F,1))'+numPoints+size(edgeMat,1);
                
                %Create faces array
                Fs=[[F(:,1)  indV_12 indV_mid indV_41];...
                    [F(:,2)  indV_23 indV_mid indV_12];...
                    [F(:,3)  indV_34 indV_mid indV_23];...
                    [F(:,4)  indV_41 indV_mid indV_34]];
                
                %Create vertex array
                Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
                
                Vm=patchCentre(F,V);
                
                Vs = [V; Vn; Vm]; %Join point sets
                
                CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1);];
                
                %% Old method
                %                 numV=size(V,1);
                %                 X=V(:,1); Y=V(:,2); Z=V(:,3);
                %
                %                 %% DERIVE EDGE INDICES AND FACE-EDGE INDEX MATRIX
                %
                %                 %Format of column index in F
                %                 EColumnInd=[(1:size(F,2)); (1:size(F,2))];
                %                 EColumnInd=[EColumnInd(2:end) EColumnInd(1)];
                %
                %                 %Derive edges matrix
                %                 E=F(:,EColumnInd)'; %Use index into F to create edges matrix
                %                 E=reshape(E,2,numel(E)/2)';
                %
                %                 E=sort(E,2); %Sort edge order
                %                 [E,~,ind2] = unique(E,'rows'); %Removing double edges, i.e. [1  4] = [4  1]
                %
                %                 Fe=reshape(1:numel(F),size(F,2),size(F,1))';
                %                 Fe=ind2(Fe);
                %
                %                 if size(F,1)==1
                %                     Fe=Fe';
                %                 end
                %
                %                 %% Calculate mid-face vertices
                %
                %                 if size(F,1)==1
                %                     XF=X(F)'; YF=Y(F)'; ZF=Z(F)';
                %                 else
                %                     XF=X(F); YF=Y(F); ZF=Z(F);
                %                 end
                %                 V_midFace=[mean(XF,2) mean(YF,2) mean(ZF,2)];
                %
                %                 %% Calculate mid-edge vertices
                %
                %                 XE=X(E); YE=Y(E); ZE=Z(E);
                %                 V_midEdge=[mean(XE,2) mean(YE,2) mean(ZE,2)];
                %
                %                 %% Create new faces matrix
                %
                %                 Vs=[V;V_midEdge;V_midFace];
                %                 startIndMidEdge=numV;
                %                 startIndMidFace=startIndMidEdge+size(V_midEdge,1);
                %                 indMidFace=startIndMidFace+1:size(Vs,1);
                %
                %                 Fs=repmat(F,4,1);
                %                 for q=1:1:4
                %                     startInd=1+(q-1)*size(F,1);
                %                     endInd=startInd-1+size(F,1);
                %
                %                     if q==1
                %                         ind4=4;
                %                     else
                %                         ind4=q-1;
                %                     end
                %
                %                     fs=[F(:,q) Fe(:,q)+numV indMidFace(:) Fe(:,ind4)+numV];
                %
                %                     Fs(startInd:endInd,:)=fs;
                %                 end
                %
                %                 %% Override input
                %                 F=Fs;
                %                 V=Vs;
            case 2
                edgeMat=[F(:,[1 2]); F(:,[3 4]);]; %Edges matrix
                Es=sort(edgeMat,2); %Sorted edges matrix
                [~,ind1,~]=unique(Es,'rows');
                edgeMat=edgeMat(ind1,:);
                
                numPoints = size(V,1);
                numEdges = size(edgeMat,1);
                
                % Get indices of the three edges associated with each face
                A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
                A = max(A,A'); %Copy symmetric
                
                %Indices for A matrix
                indA_12=F(:,1)+(F(:,2)-1)*numPoints;
                indA_34=F(:,3)+(F(:,4)-1)*numPoints;
                
                %Get indices for vertex array
                indV_12=full(A(indA_12));
                indV_34=full(A(indA_34));
                
                %Create faces array
                Fs=[[F(:,1)  indV_12 indV_34 F(:,4)];...
                    [indV_12 F(:,2)  F(:,3)  indV_34]];
                
                %Create vertex array
                Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
                
                Vs = [V; Vn; ]; %Join point sets
                
                CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1); ];
            case 3
                edgeMat=[F(:,[2 3]); F(:,[4 1]);]; %Edges matrix
                Es=sort(edgeMat,2); %Sorted edges matrix
                [~,ind1,~]=unique(Es,'rows');
                edgeMat=edgeMat(ind1,:);
                
                numPoints = size(V,1);
                numEdges = size(edgeMat,1);
                
                % Get indices of the three edges associated with each face
                A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
                A = max(A,A'); %Copy symmetric
                
                %Indices for A matrix
                indA_23=F(:,2)+(F(:,3)-1)*numPoints;
                indA_41=F(:,4)+(F(:,1)-1)*numPoints;
                
                %Get indices for vertex array
                indV_23=full(A(indA_23));
                indV_41=full(A(indA_41));
                
                %Create faces array
                Fs=[[F(:,1)  F(:,2) indV_23 indV_41];...
                    [indV_41 indV_23 F(:,3)  F(:,4)]];
                
                %Create vertex array
                Vn=0.5*(V(edgeMat(:,1),:)+V(edgeMat(:,2),:)); %new mid-edge points
                
                Vs = [V; Vn; ]; %Join point sets
                CV=[0*ones(size(V,1),1); 1*ones(size(Vn,1),1); ];
        end
        
        %% Override input
        C=repmat(C,[size(Fs,1)/size(F,1),1]);
        F=Fs;
        V=Vs;
        
    end
end

if crop2Dopt
    V=V(:,[1 2]); %Crop back to 2D
end

varargout{1}=F;
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
