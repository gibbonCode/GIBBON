function [indSeed,marchCount,seedIndex]=patchMarchCountPointSeed(F,V,indStart,n)

[~,IND_V]=patchIND(F,V);
optStruct.IND_V=IND_V;

indSeed=nan(n,1);
indSeed(1)=indStart;
for q=2:1:n
    [marchCount,~]=patchMarchCount(F,V,indSeed(~isnan(indSeed)),optStruct);
    [~,indAdd]=max(marchCount);
    indSeed(q)=indAdd;
end
[marchCount,seedIndex]=patchMarchCount(F,V,indSeed(~isnan(indSeed)),optStruct);
