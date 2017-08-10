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
% 2016/12/13 Created as improvement on (future replacement of) ind2patch
%
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
numPoints=prod(sizVox);

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
    case {'v','vu'}
        [F,C]=element2patch(F,C,'hex8');
        C_ind=repmat(C_ind,6,1);
    case 'vb'        
        [F,C]=element2patch(F,C,'hex8'); 
        C_ind=repmat(C_ind,6,1);
        [indBounary]=tesBoundary(F,numPoints);        
        F=F(indBounary,:);
        C=C(indBounary,:);
        C_ind=C_ind(indBounary,:);
end

%Get unique indices (faster than unique.m)
L=false(max(F(:)),1); 
L(F(:))=1; 
indUni=find(L);

%Fix indices in faces matrix
indFix1=1:numel(indUni);
indFix2=zeros(numPoints,1);
indFix2(indUni)=indFix1;
F=indFix2(F);

if size(F,2)==1 
    F=F'; %Tranposing F if required (possible if F is a single face)
end

if strcmp(patchType,'vu')
    %Removing double FACES
    Fs=sort(F,2); %Sort so faces with same nodes have the same rows
    [~,IND_F,IND_F_2]=unique(Fs,'rows'); %integer unique operation
    F=F(IND_F,:);
    
    %Averaring colors    
    numF=size(Fs,1); numFuni=size(F,1);
    CS=C;
    CS(CS==0)=NaN; %Set real zeros to NaN to avoid "loss" in sparse array
    sharedColourMatrixSparse=sparse(IND_F_2,1:numF,CS,numFuni,numF,numF);            
    logicColourMatrixEntry=sharedColourMatrixSparse~=0;     
    F_count=full(sum(logicColourMatrixEntry,2)); %Face counts
    C=full(nansum(sharedColourMatrixSparse,2))./F_count; %Averaging color    
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

%% Gather output

varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=C_ind;

 
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
