function [varargout]=quadSurfRemoveTwoConnect(varargin)

% function [F,V,C,indFix,L,logicTwoConnect]=quadSurfRemoveTwoConnect(F,V,C)
% -----------------------------------------------------------------------
%
%
%
%
% -----------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        C=[];
    case 3
        F=varargin{1};
        V=varargin{2};
        C=varargin{3};
end

if isempty(C)
    C=ones(size(F,1),1);
end

%%
conStruct=patchConnectivity(F,V,'vf');
conVertexFace=conStruct.vertex.face;

Eb=patchBoundary(F,V); %Boundary edges
logicBoundary=false(size(V,1),1);
logicBoundary(Eb(:))=1;
logicTwoConnect=(sum(conVertexFace>0,2)<=2) & ~logicBoundary;

indVertices=find(logicTwoConnect);
Fqn=zeros(nnz(logicTwoConnect),size(F,2));
Cqn=zeros(nnz(logicTwoConnect),size(C,2));
for q=1:1:nnz(logicTwoConnect)
    indVertexNow=indVertices(q);
    indFacesNow=conVertexFace(indVertexNow,:);
    indFacesNow=indFacesNow(indFacesNow>0);
    E=patchEdges(F(indFacesNow,:));
    logicRemove=any(E==indVertexNow,2);
    E_cand=E(~logicRemove,:);
    E1=E_cand(1,:);
    Fn=[E1 E_cand(~any(ismember(E_cand,E1),2),:)];
    
    
    Fqn(q,:)=Fn;
    Cqn(q,:)=mean(C(indFacesNow,:),1);
end
indFacesRemove=conVertexFace(logicTwoConnect,:);
indFacesRemove=indFacesRemove(indFacesRemove>0);

logicFacesKeep=true(size(F,1),1);
logicFacesKeep(indFacesRemove)=0;

F=[F(logicFacesKeep,:);Fqn];
C=[C(logicFacesKeep,:);Cqn];
L=[false(nnz(logicFacesKeep),1); true(nnz(logicTwoConnect),1)];

[F,V,indFix]=patchCleanUnused(F,V);

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=C;
varargout{4}=indFix;
varargout{5}=L;
varargout{6}=logicTwoConnect;

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
