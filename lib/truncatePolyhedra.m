function [varargout]=truncatePolyhedra(varargin)

% function [FC,VC]=truncatePolyhedra(F,V,w)
% ------------------------------------------------------------------------
% This function computes the truncated form of the input solid as specified
% by the input faces (F), vertices (V), and the truncation factor (w). The
% output is the patch data (faces, vertices, and color/face labels) for the
% truncated solid. The truncation factor ranges between 0 and 1. If
% unspecified or empty the default value is:
% w=tan(pi/(size(F,2)*2))/tan(pi/(size(F,2))); which provides the uniform
% truncation. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2018/04/03 Created to aid lattice creation
%------------------------------------------------------------------------

%%
switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        w=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        w=varargin{3};
end
        
if isempty(w)
    w=tan(pi/(size(F,2)*2))/tan(pi/(size(F,2)));
end

%%
% [NF,~,NV]=patchNormal(F,V); %Face and vertex normals

%% Get the edges matrices unique and non-unique
[E,~,indUni2]=patchEdges(F,1);
[En]=patchEdges(F,0);

%% Create the new coordinate set

Vn=zeros(numel(E),size(V,2));
for q=1:1:size(V,2) %Loop over number of dimensions
    X=V(:,q); %Coordinates
   
    XE=X(E); %Coordinates at edges
    Xm=mean(XE,2); %Mean for edges
    XEm=XE-Xm(:,ones(2,1)); %Mean subtracted to prepare for contraction
    XEm=(XEm*w)+Xm(:,ones(2,1)); %Contract and add mean
    Vn(:,q)=XEm(:); %Store coordinates in array
end

%Get new edge indices
indE=reshape((1:size(F,1)*size(F,2))',size(E,1),2);
EN=[indE(indUni2,1) indE(indUni2,2)]; 

%Flip edges which are inverted
logicFlip=E(indUni2,1)~=En(:,1); 
EN(logicFlip,:)=fliplr(EN(logicFlip,:));

%% Create "central" faces array
FN=reshape(EN',size(F,2)*2,size(F,1))';
CN=(1:1:size(FN,1))';

%% Create "corner" faces array
edgeOrder=[(2:2:size(FN,2))' [(3:2:size(FN,2))'; 1]]';
edgeOrder=edgeOrder(:);
E_FN=reshape(FN(:,edgeOrder)',2,size(FN,1)*(numel(edgeOrder)/2))';

% E_FNs=sort(E_FN,2); %Sorted so [1 4] and [4 1] are seen as the same edge
% [~,indUni1_E_FN,indUni2_E_FN]=unique(E_FNs,'rows'); %Get indices for unique edges
% E_FN=E_FN(indUni1_E_FN,:);

I=E(:);
J=(1:1:numel(I))';

S=sort(sparse(I,J,J,size(V,1),size(Vn,1),numel(J)),2,'descend');
[~,j]=find(S);
FT_ind=full(S(:,1:max(j(:))));

[~,indSort]=sort(E_FN(:,1));
E_FN=E_FN(indSort,2);

FT=FT_ind;
for q=2:1:size(FT_ind,2)
    ind=FT(:,q-1);
    FT(:,q)=E_FN(ind);
end
FT=fliplr(FT);

%% Gather output 
if size(FN,2)~=size(FT,2)
    FC={FN,FT};
    CC={zeros(size(CN)),ones(size((max(CN)+(1:size(FT,1)))'))};
else
    FC=[FN;FT];
    CC=[zeros(size(CN));ones(size((max(CN)+(1:size(FT,1)))'))];
end
VC=Vn;

varargout{1}=FC;
varargout{2}=VC;
varargout{3}=CC;
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
