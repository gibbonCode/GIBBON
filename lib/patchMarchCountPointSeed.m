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
