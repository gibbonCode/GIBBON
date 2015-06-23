function [indSeed,D_map,seedIndex]=patchMarchDistMapPointSeed(F,V,indStart,n,W)

%%
numStarts=numel(indStart);
indSeed=zeros(n,1);
indSeed(1:numStarts)=indStart; 
options.W=W; 
if n>numel(indStart)
    
    [D_map,~] = patchMarchDistMap(F,V,indStart,options);
    
    D_map_now=D_map;
    
    %Wait bar definitions
    loopRange=numStarts+1:1:n;
    numSteps=numel(loopRange);
    c=1;
    hw=waitbar(c/numSteps,['patchMarchDistMapPointSeed...',num2str(round(100.*c/numSteps)),'%']);
    for q=loopRange;
        D_map = min(D_map,D_map_now);
        [~,indSeed(q)]=max(D_map);
        options.constraint_map = D_map;
        % [D_map_now,~,~] = perform_fast_marching_mesh(V',F',indSeed(q),options);
        [D_map_now,~] = patchMarchDistMap(F,V,indSeed(q),options);
        waitbar(c/numSteps,hw,['patchMarchDistMapPointSeed...',num2str(round(100.*c/numSteps)),'%']);
        c=c+1;
    end
    close(hw);
end
[D_map,seedIndex] = patchMarchDistMap(F,V,indSeed,[]);

%%

% indSeed=zeros(n,1);
% indSeed(1)=indStart;
% [D_map,~] = patchMarchDistMap(V,F,indSeed(1),[]);
% D_map_now=D_map;
% for q=2:1:n
%     D_map = min(D_map,D_map_now);
%     [~,indSeed(q)]=max(D_map);
%     options.constraint_map = D_map;
%     [D_map_now,~,~] = patchMarchDistMap(V,F,indSeed(q),options);
% end
% [D_map,~,seedIndex] = patchMarchDistMap(V,F,indSeed,[]);

%%
% 
% indSeed=zeros(numSeeds,1);
% indSeed(1)=1; 
% 
% [D_map,~,~] = perform_fast_marching_mesh(V',F',indSeed(1),[]);
% D_map_now=D_map;
% for q=2:1:numSeeds
%     D_map = min(D_map,D_map_now);
%     [~,indSeed(q)]=max(D_map);
%     options.constraint_map = D_map;
% [D_map_now,~,~] = perform_fast_marching_mesh(V',F',indSeed(q),options);    
% end
% [D_map_now,~,seedIndex] = perform_fast_marching_mesh(V',F',indSeed);