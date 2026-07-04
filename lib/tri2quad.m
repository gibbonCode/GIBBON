function [varargout]=tri2quad(varargin)

% function [Fq,Vq,Cq]=tri2quad(Ft,Vt,centerType)
% ------------------------------------------------------------------------
% This function converts the input triangulation defined by the faces Ft
% and vertices Vt to a quandrangulation defined by the output faces Fq and
% Vq. The optional output Cq contains vertex color labels for the output
% vertices whereby Cq is 0, 1, or 3 depending on whether the output vertex
% is an original node, a new mid-edge node, or a new mid-face node. 
%
% The optional additional input centerType sets the method used. The
% following methods are supported
% centerType = 1 (DEFAULT) -> The central points are the triangle centres    
% centerType = 2 -> The central points are the triangle incentres
%
% ------------------------------------------------------------------------

%%
switch nargin
    case 2
        Ft=varargin{1};
        Vt=varargin{2};
        centerType=1; 
    case 3
        Ft=varargin{1};
        Vt=varargin{2};        
        centerType=varargin{3};         
end

%Cope with 2D input
if size(Vt,2)==2
    Vt(:,3)=0; 
end
%% Mid edge sets

edgeMat=[Ft(:,[1 2]); Ft(:,[2 3]);  Ft(:,[3 1])]; 

E_sort=sort(edgeMat,2); %Sorted edges matrix
[~,ind1,~]=unique(E_sort,'rows');
edgeMat=edgeMat(ind1,:);

numPoints = size(Vt,1);
numEdges = size(edgeMat,1);

% Get indices of the four edges associated with each face
A = sparse(edgeMat(:,1),edgeMat(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
A = max(A,A'); %Copy symmetric

%Indices for A matrix
indA_12=Ft(:,1)+(Ft(:,2)-1)*numPoints;
indA_23=Ft(:,2)+(Ft(:,3)-1)*numPoints;
indA_31=Ft(:,3)+(Ft(:,1)-1)*numPoints;

%Get indices for vertex array
indV_12=full(A(indA_12));
indV_23=full(A(indA_23));
indV_31=full(A(indA_31));

%% Mid face

indV_midFace=(1:1:size(Ft,1))';
indOffset=numPoints+size(edgeMat,1);

indV_midFace123=(indV_midFace((1-1)*size(Ft,1)+(1:size(Ft,1))))+indOffset;

%% Create quad faces array

 Fq=[Ft(:,1) indV_12 indV_midFace123 indV_31;... %Corner quad 1
    indV_12  Ft(:,2) indV_23 indV_midFace123;... %Corner quad 2
     indV_midFace123 indV_23 Ft(:,3) indV_31;... %Corner quad 3     
    ];

%% Create vertex arrays

%new mid-edge points
Vn=0.5*(Vt(edgeMat(:,1),:)+Vt(edgeMat(:,2),:));
     
switch centerType
    case 1 %Use mean of triangle coordinates for central point                       
        Vm=zeros(size(Ft,1),3);
        for q=1:1:size(Vt,2)
            X=Vt(:,q);
            if size(Ft,1)==1
                Vm(:,q)=mean(X(Ft)',2);
            else
                Vm(:,q)=mean(X(Ft),2);
            end
        end        
    case 2 %Use triangle incentres for central point                
        Vm=triIncenter(Ft,Vt);
end

Vq=[Vt; Vn; Vm]; %Join point sets

CVq=[0*ones(size(Vt,1),1); 1*ones(size(Vn,1),1); 2*ones(size(Vm,1),1);]; %Vertex colors for original, mid-edge, and mid-face nodes

%% Collect ouput

varargout{1}=Fq;
varargout{2}=Vq;
varargout{3}=CVq;

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
