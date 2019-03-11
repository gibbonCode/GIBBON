function [varargout]=subEdge(varargin)

%%

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        n=1;
        uniqueOpt=1;
    case 3
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        uniqueOpt=1;
    case 4
        F=varargin{1};
        V=varargin{2};
        n=varargin{3};
        uniqueOpt=varargin{4};
end

%%

%Check if n can be achieved through splitting
nSplitIterations=log2(n+1); %Check for integer solution
logicInteger=abs(round(nSplitIterations)-nSplitIterations)<eps(nSplitIterations);

if uniqueOpt==1 && logicInteger
    subMethod='split';
else
    subMethod='seed';
end

switch subMethod
    case 'split' %iteratively split edges (no double points created, unique operation avoided)
        Fs=F; Vs=V;
        for q=1:1:nSplitIterations
            F=Fs; V=Vs;
            
            Fs=zeros(size(Fs,1),size(Fs,2)*2);
            Fs(:,1:2:end)=F;
            
            [E]=patchEdges(F,1);
            
            numPoints = size(V,1);
            numEdges = size(E,1);
            
            % Get indices of the edges associated with each face
            A = sparse(E(:,1),E(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);
            A = max(A,A'); %Copy symmetric
            
            indFs=2:2:size(Fs,2);
            for qp=1:1:size(F,2)
                if qp<size(F,2)
                    indA_now=F(:,qp)+(F(:,qp+1)-1)*numPoints;
                else
                    indA_now=F(:,qp)+(F(:,1)-1)*numPoints;
                end
                
                indV_now=full(A(indA_now));
                Fs(:,indFs(qp))=indV_now;
            end
            
            %Create vertex array
            Vn=0.5*(V(E(:,1),:)+V(E(:,2),:)); %new mid-edge points
            Vs = [V; Vn]; %Join point sets
        end
    case 'seed' %Seed edge points and remove doubles (more memory intensive)
        %Edges matrix
        [E]=patchEdges(F,0);
        
        Vs=zeros(size(E,1)*(n+1),size(V,2));
        for q=1:1:size(V,2)
            X=V(:,q);
            XE=X(E);
            X_add=linspacen(XE(:,1),XE(:,2),n+2);
            X_add=X_add(:,1:end-1)';
            Vs(:,q)=X_add(:);
        end
        
        ind=(1:1:size(Vs,1))';
        Fs=reshape(ind,[(n+1)*size(F,2),size(Vs,1)./((n+1)*size(F,2))])';
        
        %Merge non-unique nodes
        [Fs,Vs]=mergeVertices(Fs,Vs);
        
end
 
%% Collect output
varargout{1}=Fs;
varargout{2}=Vs;


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
