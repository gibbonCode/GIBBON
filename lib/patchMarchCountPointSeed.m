function [indSeed,marchCount,seedIndex]=patchMarchCountPointSeed(F,V,indStart,n)

[~,IND_V]=tesIND(F,V,0);
optStruct.IND_V=IND_V;

numStarts=numel(indStart);
indSeed=nan(n,1);
indSeed(1:numStarts)=indStart; 

for q=numStarts:1:n+1
    [marchCount,seedIndex]=patchMarchCount(F,V,indSeed(~isnan(indSeed)),optStruct);
    [~,indAdd]=max(marchCount);
    if q<=n
        indSeed(q)=indAdd;
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
