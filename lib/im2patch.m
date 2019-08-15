function [varargout]=im2patch(varargin)

% function [F,V,C,C_ind]=im2patch(M,indPatch,patchType,v)
% ------------------------------------------------------------------------
%
% This function generates patch data (faces 'F', vertices 'V' and color
% data 'C') for 2D or 3D images. The patches are only generated for the
% voxels specified by the logic or linear indices in 'indPatch'. The
% variable 'patchType' indicates the type of patch: 
%
% 'v'              Voxel patch data whereby each voxel has 8 potentially
%                  shared vertices and 6 unshared faces. 
% 'vu'             Voxel patch data whereby each voxel has 8 potentially
%                  shared vertices and 6 potentially shared faces. Colors
%                  for shared faces is averaged. 
% 'vb'             Voxel patch data whereby only unshared and boundary
%                  faces are displayed. This is the most lightweight voxel
%                  patch type. 
% 'si', 'sj', 'sk' Mid-voxel slice patch data for the i, j and k direction
%                  respectively. Vertices are shared. 
% 'h'              Creates a hexahedral element description of the voxels
%                  (e.g. a 8-column hexahedral element matrix). 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2016/12/13 Created as improvement on (and future replacement of) ind2patch
% 2018/12/19 Added patchCleanUnused to remove unused vertices
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        M=varargin{1};
        indPatch=true(size(M));
        patchType=[];
        v=[];
    case 2
        M=varargin{1};
        indPatch=varargin{2};
        patchType=[];
        v=[];
    case 3
        M=varargin{1};
        indPatch=varargin{2};
        patchType=varargin{3};
        v=[];
    case 4
        M=varargin{1};
        indPatch=varargin{2};
        patchType=varargin{3};
        v=varargin{4};
    otherwise
        error('Wrong number of input arguments')
end

if isempty(patchType)
    patchType='vb';
end

if isempty(indPatch)
    varargout{1}=[];
    varargout{2}=[];
    varargout{3}=[];
    varargout{4}=[];
    return
end

%Convert logical to linear indices
if islogical(indPatch) 
    indPatch=find(indPatch); 
end
indPatch=indPatch(:);

%Deal with complex data
if ~isreal(M)
    M=abs(M);
    disp('Warning: Complex data, magnitudes of complex numbers are used for colors');
end

%%
siz=size(M);
if numel(siz)==2
    siz(3)=1;
end

sizVox=siz+1;

[I,J,K]=ind2sub(siz,indPatch);
% numPoints=prod(sizVox);

%Prepare colors
C=double(M(indPatch));
C_ind=indPatch;

%Build faces matrix
switch patchType
    case {'h','v','vb','vu'}
        EI=[I I+1 I+1 I I I+1 I+1 I];
        EJ=[J J J+1 J+1 J J J+1 J+1];
        EK=[K K K K K+1 K+1 K+1 K+1];                       
    case 'si'
        EI=[I  I   I   I];
        EJ=[J  J+1 J+1 J];
        EK=[K  K   K+1 K+1];
    case 'sj'
        EI=[I  I+1 I+1 I];
        EJ=[J  J   J   J];
        EK=[K  K   K+1 K+1];
    case 'sk'
        EI=[I  I+1 I+1 I];
        EJ=[J  J   J+1 J+1];
        EK=[K  K   K   K ];  
    otherwise
        error('invalid patch type');
end
F=sub2ind(sizVox,EI,EJ,EK);                        

switch patchType
    case 'h'
        F=F(:,[4:-1:1 8:-1:5]); %Invert voxel face order to form elements
end

%Alter face description for voxels
switch patchType
    case {'v','vu','vb'}
        [F,C]=element2patch(F,C,'hex8');
        C_ind=repmat(C_ind,6,1);        
end

%Get unique indices (faster than unique.m)
L=false(max(F(:)),1); 
L(F(:))=1; 
indUni=find(L);

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(max(F(:)),1);
indFix2(indUni)=indFix1;
F=indFix2(F);

if size(F,2)==1 
    F=F'; %Tranposing F if required (possible if F is a single face)
end

switch patchType
    case {'vb'} %Get only boundary faces
        [indBounary]=tesBoundary(F,max(F(:)));
        F=F(indBounary,:);
        C=C(indBounary,:);
        C_ind=C_ind(indBounary,:);
    case {'vu'} %Removing double FACES
        Fs=sort(F,2); %Sort so faces with same nodes have the same rows
        [~,IND_F,IND_F_2,F_count]=cunique(Fs,'rows'); %get indices for unique faces   
        F_count=F_count(IND_F);
        
        F=F(IND_F,:);
        
        %Averaring colors
        numF=size(Fs,1); numFuni=size(F,1);        
        C(isnan(C))=0; % Replace NaN's by zeros
        sharedColourMatrixSparse=sparse(IND_F_2,1:numF,C,numFuni,numF,numF);        
        C=full(sum(sharedColourMatrixSparse,2))./F_count; %Averaging color
        
        C_ind=sort(sparse(IND_F_2,1:numF,C_ind,numFuni,numF,numF),2);
        C_ind=full(C_ind(:,end-max(F_count)+1:end));
end

%Create coordinate set
[Iv,Jv,Kv]=ind2sub(sizVox,indUni);
Iv=Iv-0.5; Jv=Jv-0.5; Kv=Kv-0.5;

%Shift slice types to middle of voxel
switch patchType        
    case 'si'
        Iv=Iv+0.5;         
    case 'sj'
        Jv=Jv+0.5;         
    case 'sk'        
        Kv=Kv+0.5; 
end

if ~isempty(v)
    [Xv,Yv,Zv]=im2cart(Iv,Jv,Kv,v);
    V=[Xv(:) Yv(:) Zv(:)];
else
    V=[Jv Iv Kv]; %N.B. I and J direction are switched
end

[F,V]=patchCleanUnused(F,V); % Remove unused vertices

%% Gather output

varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=C_ind;
 
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
