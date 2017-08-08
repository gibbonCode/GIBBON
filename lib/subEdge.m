function [Fs,Vs]=subEdge(varargin)

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
        numDigitKeep=5;
        [~,ind1,ind2]=unique(pround(Vs,numDigitKeep),'rows');
        Vs=Vs(ind1,:);
        
        %Treat special case of 1 face
        if size(Fs,1)>1
            Fs=ind2(Fs);
        else
            Fs=ind2(Fs)';
        end
end
 
%% 
% ********** _license boilerplate_ **********
% 
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
