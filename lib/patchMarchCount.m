function [marchCount,seedIndex]=patchMarchCount(F,V,indStart,optStruct)

if ~isempty(optStruct)
    if isfield(optStruct,'IND_V')
        IND_V=optStruct.IND_V;
    else
        [~,IND_V,~]=tesIND(F,V,0);
    end
else
    [~,IND_V,~]=tesIND(F,V,0);
end

marchCount=zeros(size(V,1),1);
indGroupNew=indStart;

seedIndex=zeros(size(V,1),1);
seedIndex(indStart)=indStart;
seedIndex_IND_V=IND_V;
seedIndex_IND_V(IND_V>0)=seedIndex(IND_V(IND_V>0));
seedIndex_IND_V(indStart,1)=indStart;
seedIndex_IND_V=max(seedIndex_IND_V,[],2);
seedIndex_IND_V=seedIndex_IND_V(:,ones(1,size(IND_V,2)));

logicGroup=false(size(V,1),1);
logicGroup(indStart)=1;
logicGroupNew=false(size(V,1),1);
logicGroupNew(indGroupNew)=1;
countLevel=1;
numPoints=size(V,1);

indAll=1:1:size(V,1);

while 1
    marchCount(logicGroupNew)=countLevel;
    
    indUsed=indAll(logicGroupNew); %KMM edit 2014/05/26 to ensure points are used once
    logicUsed=ismember(IND_V,indUsed); %KMM edit 2014/05/26 to ensure points are used once
    IND_V(logicUsed)=0; %KMM edit 2014/05/26 to ensure points are used once
    
    %Recompose current group
    indGroupNew=IND_V(logicGroupNew,:); %Current and new group members and zeros
    Lv=indGroupNew>0;
    indGroupNew=indGroupNew(Lv); %Group members without zeros
    
    seedIndexNew=seedIndex_IND_V(logicGroupNew,:);
    seedIndex_IND_V(indGroupNew,1)=seedIndexNew(Lv);
    seedIndex=max(seedIndex_IND_V,[],2);
    seedIndex_IND_V=seedIndex(:,ones(1,size(IND_V,2)));
    
    logicGroupNew=false(numPoints,1);
    logicGroupNew(indGroupNew)=1;
    logicGroupNew=logicGroupNew & ~logicGroup; %Remove already existing
    
    logicGroup(indGroupNew)=1;
    countLevel=countLevel+1;
    
    if nnz(marchCount)==numPoints
        break
    end
end

 
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
