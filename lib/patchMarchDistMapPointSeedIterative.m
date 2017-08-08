function [indSeed,D_map,seedIndex]=patchMarchDistMapPointSeedIterative(F,V,indStart,numSeeds,W,optStruct)

%%

if isempty(W)
    W=ones(size(V,1),1);
end

if isempty(indStart)
    indStart=1;
end

if ~isfield(optStruct,'numIterationsMax')
    optStruct.numIterationsMax=200;
end

if ~isfield(optStruct,'toleranceLevel')
    optStruct.toleranceLevel=1e-4;
end


if isfield(optStruct,'waitBarOn')
    waitBarOn=optStruct.waitBarOn;
else
    waitBarOn=0;
end


%%

if numSeeds<numel(indStart)
    error('Number of seeds should be equal to or larger than the number of start points')
end

indSeed=indStart;

if waitBarOn
    hw = waitbar(0,'Please wait...');
end

for q=numel(indStart):1:numSeeds
    if waitBarOn
        waitbar(q/numSeeds,hw,['patchMarchDistMapPointSeedIterative ',num2str(round(q/numSeeds*100)),'% complete']);
    end
    
    [D_map] = patchMarchDistMapIterative(F,V,indSeed(1:q),W,optStruct); %New distance map

    if q>=numel(indStart) && q<numSeeds
        [~,indSeed(q+1)]=max(D_map);
    end    
end

if waitBarOn
    close(hw);
end

[~,~,seedIndex]=patchMarchCountPointSeed(F,V,indSeed,numel(indSeed));


% %%
% 
% if numSeeds<numel(indStart)
%     error('Number of seeds should be equal to or larger than the number of start points')
% end
% 
% [E]=patchEdgeLengths(F,V);
% eMin=min(E(:));
% 
% seedIndex=zeros(size(V,1),1);
% indSeed=indStart;
% 
% if waitBarOn
%     hw = waitbar(0,'Please wait...');
% end
% 
% for q=1:1:numSeeds
%     if waitBarOn
%         waitbar(q/numSeeds,hw,['patchMarchDistMapPointSeedIterative ',num2str(round(q/numSeeds*100)),'% complete']);
%     end
%     [D_map] = patchMarchDistMapIterative(F,V,indSeed(1:q),W,optStruct); %New distance map
%     if q==1
%         seedIndex=ones(size(V,1),1)*indSeed(q);
%     else
%         d=(D_map-D_map_now);
%         logicSeedIndex=d<-(eMin/100);
%         seedIndex(logicSeedIndex)=indSeed(q);
%     end
%     if q>=numel(indStart) && q<numSeeds
%         [~,indSeed(q+1)]=max(D_map);
%     end
%     D_map_now=D_map;
% end
% 
% if waitBarOn
%     close(hw);
% end
 
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
