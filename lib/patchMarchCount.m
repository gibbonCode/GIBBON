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

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2017  Kevin Mattheus Moerman
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
