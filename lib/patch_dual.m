function [varargout]=patch_dual(varargin)

% function [Vd,Fd,Fds]=patch_dual(V,F,fixBoundaryOption)
% ------------------------------------------------------------------------
% Computes the dual of the input tesselation defined by the vertices V and
% faces F. 
%
%
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2013/04/28
% 2015/07/07 Updated for GIBBON
% 2018/04/06 Added color output
% 2018/04/06 Improved handling of edge cells
%------------------------------------------------------------------------
%%

switch nargin
    case 2
        V=varargin{1};
        F=varargin{2};
        fixBoundaryOption=1;
    case 3
        V=varargin{1};
        F=varargin{2};
        fixBoundaryOption=varargin{3};
end

%Cope with 2D input
if size(V,2)==2
    V(:,3)=0;
end

%%

if iscell(F)
%     numFacesperSet=cellfun(@(x) size(x,1),F);
    
    N=[];
    vertexFaceConnectivity=[];
    numFacesPrevous=0;
    for q=1:1:numel(F)
        %The current face set
        f=F{q};
        
        %Get indices of face neighbourhood
        [vertexFaceConnectivity_f]=tesIND(f,V,0);    
        
        %Offset face indices by previous 
        vertexFaceConnectivity_f(vertexFaceConnectivity_f>0)=vertexFaceConnectivity_f(vertexFaceConnectivity_f>0)+numFacesPrevous;
        
        %Get face-vertex connectivity
        vertexFaceConnectivity=[vertexFaceConnectivity vertexFaceConnectivity_f];
        
        %Derive vertex normals
        N=[N; patchNormal(f,V)];

        numFacesPrevous=numFacesPrevous+size(f,1);        
    end
    vertexFaceConnectivity=sort(vertexFaceConnectivity,2);
    [~,J,~] = find(vertexFaceConnectivity);
    vertexFaceConnectivity=vertexFaceConnectivity(:,min(J):end);
    vertexFaceConnectivity=full(vertexFaceConnectivity);
    
    L=vertexFaceConnectivity>0;
    
    Nv=ones(size(V,1),size(N,2));
    for q=1:1:size(N,2)
        nf=N(:,q);
        nv=nan(size(vertexFaceConnectivity));
        nv(L)=nf(vertexFaceConnectivity(L));
        nv=gnanmean(nv,2);
        Nv(:,q)=nv;
    end
else
    %Get indices of face neighbourhood
    [vertexFaceConnectivity]=tesIND(F,V,0);
    
    %Get input tesselation vertex-normals
    [~,~,Nv]=patchNormal(F,V);
end

%Create central face coordinates
if iscell(F)
    VF=[];
    for q=1:1:numel(F)        
        VF=[VF; patchCentre(F{q},V)];
    end
else
    VF=patchCentre(F,V);
end

if fixBoundaryOption>0
    
    [Eb]=patchBoundary(F,V);
    if ~isempty(Eb) 
        switch fixBoundaryOption
            case 1 %Include original boundary vertices
                ind_Eb=unique(Eb(:));
                
                %Create mid-edge vertices
                VE=zeros(size(Eb,1),size(V,2));
                for q=1:1:size(V,2)
                    X=V(:,q);
                    VE(:,q)=mean(X(Eb),2);
                end
                
                Vd=[VF; VE; ];
                
                %Add new points in IND_F
                vertexFaceConnectivity(Eb(:,1),end+1)=(1:size(VE,1))'+size(VF,1);
                ind1=vertexFaceConnectivity(ind_Eb,end);
                vertexFaceConnectivity(Eb(:,2),end+1)=(1:size(VE,1))'+size(VF,1);
                ind2=vertexFaceConnectivity(ind_Eb,end);
                vertexFaceConnectivity(ind_Eb,end+1)=(1:1:numel(ind_Eb))'+size(VF,1)+size(VE,1);
                
                VB=(Vd(ind1,:)+Vd(ind2,:))/2; %Use this average for now to avoid inverted elements VB=V(ind_Eb,:);
                                
                %Collect point sets
                Vd=[VF; VE; VB];
                
                vertexFaceConnectivity=sort(vertexFaceConnectivity,2,'descend');
                indMax=find(sum(vertexFaceConnectivity,1)==0,1);
                vertexFaceConnectivity=vertexFaceConnectivity(:,1:indMax);
                
                replaceVb=1; %Replace on so average is replaced later by real boundary vertices
            case 2 % Do not include original boundary vertices
                %Create mid-edge vertices
                VE=zeros(size(Eb,1),size(V,2));
                for q=1:1:size(V,2)
                    X=V(:,q);
                    VE(:,q)=mean(X(Eb),2);
                end
                
                %Collect point sets
                Vd=[VF; VE; ];
                
                %Add new points in IND_F
                vertexFaceConnectivity(Eb(:,1),end+1)=(1:size(VE,1))'+size(VF,1);
                vertexFaceConnectivity(Eb(:,2),end+1)=(1:size(VE,1))'+size(VF,1);                
                
                vertexFaceConnectivity=sort(vertexFaceConnectivity,2,'descend');
                indMax=find(sum(vertexFaceConnectivity,1)==0,1);
                vertexFaceConnectivity=vertexFaceConnectivity(:,1:indMax);
                
                replaceVb=0;
        end
%         %Create boundary vertices
%         ind_Eb=unique(Eb(:));
%         VB=V(ind_Eb,:);
%         
%         %Create mid-edge vertices
%         VE=zeros(size(Eb,1),size(V,2));
%         for q=1:1:size(V,2)
%             X=V(:,q);
%             VE(:,q)=mean(X(Eb),2);
%         end
%         
%         %Collect point sets
%         Vd=[VF; VE; VB];
%         
%         %Add new points in IND_F
%         IND_F(Eb(:,1),end+1)=(1:size(VE,1))'+size(VF,1);
%         IND_F(Eb(:,2),end+1)=(1:size(VE,1))'+size(VF,1);
%         IND_F(ind_Eb,end+1)=(1:1:numel(ind_Eb))'+size(VF,1)+size(VE,1);
%         
%         IND_F=sort(IND_F,2,'descend');
%         indMax=find(sum(IND_F,1)==0,1);
%         IND_F=IND_F(:,1:indMax);

    else
        Vd=VF;
        replaceVb=0;
    end
else
    Vd=VF;
    replaceVb=0;
end

%Creating arrays for faces
[I,J,v] = find(vertexFaceConnectivity);

Xfd=accumarray({I,J},Vd(v,1),size(vertexFaceConnectivity),[],NaN); 
Yfd=accumarray({I,J},Vd(v,2),size(vertexFaceConnectivity),[],NaN); 
Zfd=accumarray({I,J},Vd(v,3),size(vertexFaceConnectivity),[],NaN); 

Xfd_mean=gnanmean(Xfd,2); 
Yfd_mean=gnanmean(Yfd,2); 
Zfd_mean=gnanmean(Zfd,2); 

Xfd=Xfd-Xfd_mean(:,ones(1,(size(Xfd,2))));
Yfd=Yfd-Yfd_mean(:,ones(1,(size(Yfd,2))));
Zfd=Zfd-Zfd_mean(:,ones(1,(size(Zfd,2))));

%Determine face order
Xfdn=Xfd; Yfdn=Yfd; Zfdn=Zfd;
for q=1:1:size(Xfd,1)
    Vc=[Xfd(q,:); Yfd(q,:); Zfd(q,:)]';    
    logicValid=~any(isnan(Vc),2);
    [R]=pointSetPrincipalDir(Vc(logicValid,:)); %Fit local coordinate system with 3rd direction pointing outward of local planar-ish region
    Vcn=Vc*R; %Rotate coordinate set to prepare for 2D Delaunay based triangulation    
    Xfdn(q,:)=Vcn(:,1);
    Yfdn(q,:)=Vcn(:,2);
    Zfdn(q,:)=Vcn(:,3);    
end

%Fix face order
theta_n=atan2(Yfdn,Xfdn); 
[~,J_sort]=sort(theta_n,2);
I_sort=(1:1:size(J_sort,1))'; 
I_sort=I_sort(:,ones(1,size(J_sort,2)));
IND_sort = sub2ind(size(J_sort),I_sort,J_sort);

% Creating faces matrix
Fds=vertexFaceConnectivity(IND_sort);

%Splitting up into seperate face types and fix face normals
n_sum=sum(Fds>0,2);
faceTypesNum=unique(n_sum);
faceTypesNum=faceTypesNum(faceTypesNum>0);
Fd=cell(1,numel(faceTypesNum));
C=(1:1:size(V,1))';
Cd=cell(1,numel(faceTypesNum));
for q=1:1:numel(faceTypesNum)
    %Get faces
    logicFacesNow=(n_sum==faceTypesNum(q)); %logic for current faces    
    F_now=Fds(logicFacesNow,1:faceTypesNum(q)); %The current face set
    
    %Flip face orientation if required
    Nv_now=Nv(logicFacesNow,:); %Appropriate face normals based on input mesh    
    [N_now]=patchNormal(F_now,Vd); %Current face normals        
    logicFlip=dot(Nv_now,N_now,2)<0; %dot product is negative if face is inverted
    F_now(logicFlip,:)=fliplr(F_now(logicFlip,:)); %Flip faces
  
    indFaces=find(logicFacesNow);
    indFlip=indFaces(logicFlip);
    Fds(indFlip,1:faceTypesNum(q))=fliplr(Fds(indFlip,1:faceTypesNum(q)));
    
    Fd{q}=F_now; %Store faces in cell array    
    Cd{q}=C(logicFacesNow);
end

if replaceVb==1 && fixBoundaryOption~=2
    Vd(end-numel(ind_Eb)+1:end,:)=V(ind_Eb,:);
end

varargout{1}=Vd; %Vertices
varargout{2}=Fd; %Cell array of faces
varargout{3}=Fds; %Array of face indices (with zeros)
varargout{4}=Cd; %Colors
 
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
