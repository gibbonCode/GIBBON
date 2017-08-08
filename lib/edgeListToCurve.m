function [indList]=edgeListToCurve(E)

if size(E,1)==1
    indList=E;
else
    [~,indV,~]=tesIND(E,[],0);
    logicEndPoints=sum(indV>0,2)==1;
    
    if nnz(logicEndPoints)==0 %closed loop
        indStartPoint=E(1,1);
    else
        indStartPoint=find(logicEndPoints,1);
    end
    indStartEdge=find(any(E==indStartPoint,2),1);
    E_now=E(indStartEdge,:);
    
    indUni=unique(E(:));
    
    ind1=E(indStartEdge,E_now==indStartPoint);
    ind2=E(indStartEdge,E_now~=indStartPoint);
    indList=ones(1,numel(indUni));
    indList(1)=ind1;
    indList(2)=ind2;
    q=3;
    Es=E;
    Es(indStartEdge,:)=NaN;
    while 1
        [indE_next,~]=find(Es==ind2,1);
        E_now=Es(indE_next,:);
        Es(indE_next,:)=NaN;
        ind3=E_now(E_now~=ind2);
        
        
        indList(q)=ind3;
        ind2=ind3;
        if nnz(isnan(Es))==numel(Es)
            break
        end
        q=q+1;
    end
    
    %% Fix order
    logicFlip=(E(:,1)==indList(1)) & (E(:,2)==indList(2));
    if ~any(logicFlip)
        indList=flipud(indList); %invert curve to conform to edge directions
    end
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
