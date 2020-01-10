function [varargout]=triSurfRemoveThreeConnect(varargin)

% function [Ft,Vt,Ct,L]=triSurfRemoveThreeConnect(Fd,Vd,Cd)
% ------------------------------------------------------------------------
% In a surface triangulation "3-connected" locations often contain poor
% quality triangles of a locally smaller area then the rest of the surface.
% Smoothening does not resolve this issue since the quality is not great
% improved even after vertex is at the centre of its neighbouring nodes.
% Hence the function triSurfRemoveThreeConnect instead removes the central
% nodes and groups the affected triangles into a single triangle.
% The input sets Fd, Vd and Cd represent the input faces, vertices and face
% colours respectively. The output sets Ft, Vt and Ct represent the fixed
% faces, vertices and face colours respectively. In addition a logic L can
% be output which defines the affected vertices in Vd. The output may
% consist of [Ft,Vt] or [Ft,Vt,Ct] or [Ft,Vt,Ct,L].
% The last nnz(L)/3 faces in Ft represent the fixed faces.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/06/03 Created
% 2015/05/01 Updated with varagin type input
%------------------------------------------------------------------------

%% Parse input
Fd=varargin{1};
Vd=varargin{2};

if nargin==3
    Cd=varargin{3};
else
    Cd=[];
end

if isempty(Cd)
    Cd=(1:1:size(Fd,1))';
end

%%

%Get patch face/vertex connectivity matrices
[IND_F,IND_V]=patchIND(Fd,Vd,1);

%Count point connectivity
numFriends=sum(IND_V>0,2); %Number of vertex neighbours
logicThree=numFriends==3; %Logic for vertices with only three connected neighbours

%Remove boundary points from list
TR=triangulation(Fd,Vd);
indFree=freeBoundary(TR);
indFree=unique(indFree(:));
if ~isempty(indFree)
    logicThree(indFree)=0;
end

if nnz(logicThree)>0
    
    %IND_V subset
    IND_V_three=IND_V(logicThree,:);
    
    %Snap vertices to mean of neighbours such that vertex normal should
    %coincide with new surface normal
    logicValid=IND_V_three>0;
    Xp=NaN(size(IND_V_three));
    Yp=NaN(size(IND_V_three));
    Zp=NaN(size(IND_V_three));
    Xp(logicValid)=Vd(IND_V_three(logicValid),1);
    Yp(logicValid)=Vd(IND_V_three(logicValid),2);
    Zp(logicValid)=Vd(IND_V_three(logicValid),3);
    Vp=[gnanmean(Xp,2) gnanmean(Yp,2) gnanmean(Zp,2)];
    Vd(logicThree,:)=Vp; %Replace points
    
    %Get new faces
    IND_F_three=IND_F(logicThree,:);
    indFacesThree=IND_F_three(IND_F_three>0);
    logicFacesThree=false(size(Fd,1),1);
    logicFacesThree(indFacesThree)=1;
    IND_V_three=sort(IND_V_three,2);
    if size(IND_V_three,2)>3
        IND_V_three=IND_V_three(:,end-2:end);
    end
    F_new=IND_V_three;
    
    %Fix face orientation based on normals
    indUsed=unique(F_new(:));
    V_new=Vd(indUsed,:);
    indFix=zeros(size(Vd,1),1);
    indFix(indUsed)=1:numel(indUsed);
    F_new_fix=indFix(F_new);
    if size(F_new,1)==1
        F_new_fix=F_new_fix';
    end    
    [~,~,Nd]=patchNormal(Fd,Vd);
    [N_new]=patchNormal(F_new_fix,V_new);
    N_sum_mag=sqrt(sum((Nd(logicThree,:)+N_new).^2,2));
    logicFlip=N_sum_mag<1;
    F_new(logicFlip,:)=fliplr(F_new(logicFlip,:));
    
    %Get color information for new faces
    logicValid=IND_F_three>0;
    
    C_F_three=nan(size(IND_F_three));
    C_F_three(logicValid)=Cd(IND_F_three(logicValid));
    C_new=gnanmean(C_F_three,2);
    
    %Faces to keep
    F_keep=Fd(~logicFacesThree,:);
    C_keep=Cd(~logicFacesThree,:);
    
    %Join and create vertex/face/color sets
    Ft=[F_keep;F_new];
    Ct=[C_keep;C_new];
    
    indUsed=unique(Ft(:));
    Vt=Vd(indUsed,:);
    indFix=zeros(size(Vd,1),1);
    indFix(indUsed)=1:numel(indUsed);
    Ft=indFix(Ft);
    
else
    Ft=Fd;
    Vt=Vd;
    Ct=Cd;
    logicFacesThree=false(size(Ft,1),1);
end

switch nargout
    case 2
        varargout{1}=Ft;
        varargout{2}=Vt;
    case 3
        varargout{1}=Ft;
        varargout{2}=Vt;
        varargout{3}=Ct;
    case 4
        varargout{1}=Ft;
        varargout{2}=Vt;
        varargout{3}=Ct;
        varargout{4}=logicFacesThree;
    otherwise
        error('Wrong number of output arguments');
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
