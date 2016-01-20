function [D_map,seedIndex]=patchMarchDistMap(F,V,indStart,options)

% try    
    [D_map,~,seedIndex] = perform_fast_marching_mesh(V',F',indStart,options); % See: www.numerical-tours.com    
% catch
%     warning('Missing external functions using experimental code instead!')
%     
%     numPoints=size(V,1);
%     
%     [~,IND_V]=patchIND(F,V);
%     
%     marchCount=zeros(numPoints,1);
%     indGroupNew=indStart;
%     
%     seedIndex=zeros(numPoints,1);
%     seedIndex(indStart)=indStart;
%     seedIndex_IND_V=IND_V;
%     seedIndex_IND_V(IND_V>0)=seedIndex(IND_V(IND_V>0));
%     seedIndex_IND_V(indStart,1)=indStart;
%     seedIndex_IND_V=max(seedIndex_IND_V,[],2);
%     seedIndex_IND_V=seedIndex_IND_V(:,ones(1,size(IND_V,2)));
%     
%     logicGroup=false(numPoints,1);
%     logicGroup(indStart)=1;
%     indGroup=find(logicGroup);
%     logicGroupNew=false(numPoints,1);
%     logicGroupNew(indGroupNew)=1;
%     countLevel=1;
%     
%     %Coords of vertices
%     X=V(:,1);
%     Y=V(:,2);
%     Z=V(:,3);
%     
%     %Coords of neighbouring vertices in IND_V "form"
%     IND_V_X=nan(size(IND_V));
%     IND_V_X(IND_V>0)=X(IND_V(IND_V>0));
%     
%     IND_V_Y=nan(size(IND_V));
%     IND_V_Y(IND_V>0)=Y(IND_V(IND_V>0));
%     
%     IND_V_Z=nan(size(IND_V));
%     IND_V_Z(IND_V>0)=Z(IND_V(IND_V>0));
%     
%     %Coords of vertices expanded to IND_V "form"
%     X=X(:,ones(size(IND_V,2),1));
%     Y=Y(:,ones(size(IND_V,2),1));
%     Z=Z(:,ones(size(IND_V,2),1));
%     
%     D_map=zeros(numPoints,1);
%     D_map(indStart)=0;
%     
%     
%     while 1
%         
%         marchCount(logicGroupNew)=countLevel;
%         
%         %Compose the current group
%         indGroupNew=IND_V(logicGroupNew,:); %Previous and and new group members and zeros
%         logicVisited=ismember(indGroupNew,indGroup);
%         indGroupNew(logicVisited)=0; %New members not visited yet
%         Lv=indGroupNew>0; %Logic for valid entries
%         indGroupNew=indGroupNew(Lv); %Group members without zeros
%         indGroupNew=indGroupNew(:); %Place on column
%         
%         %Derive distance metrics
%         DX=IND_V_X(logicGroupNew,:)-X(logicGroupNew,:);
%         DY=IND_V_Y(logicGroupNew,:)-Y(logicGroupNew,:);
%         DZ=IND_V_Z(logicGroupNew,:)-Z(logicGroupNew,:);
%         D_sub=sqrt(DX.^2+DY.^2+DZ.^2);
%         D_sub(~Lv)=NaN;
%         D_prev=D_map(:,ones(size(IND_V,2),1));
%         D_prev=D_prev(logicGroupNew,:);
%         D_sum=D_sub+D_prev;
%         D_map(indGroupNew)=D_sum(~isnan(D_sum));
%         
%         %SeedIndex
%         seedIndexNew=seedIndex_IND_V(logicGroupNew,:);
%         seedIndex_IND_V(indGroupNew,1)=seedIndexNew(Lv);
%         seedIndex=max(seedIndex_IND_V,[],2);
%         seedIndex_IND_V=seedIndex(:,ones(1,size(IND_V,2)));
%         
%         logicGroupNew=false(numPoints,1);
%         logicGroupNew(indGroupNew)=1;
%         logicGroupNew=logicGroupNew & ~logicGroup; %Remove already existing
%         
%         logicGroup(indGroupNew)=1;
%         indGroup=find(logicGroup);
%         countLevel=countLevel+1;
%         
%         if nnz(marchCount)==numPoints
%             break
%         end
%     end
% end

