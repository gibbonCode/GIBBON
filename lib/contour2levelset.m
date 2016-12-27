function [K]=contour2levelset(M,v,Vcs)


%%
resampleOpt=0;
siz=size(M); %Size of the image dataset

%Get image coordinates
[J,I]=meshgrid(1:1:siz(2),1:1:siz(1));

%Convert to cartesian coordinates using voxel size if provided
[x,y,~]=im2cart(I,J,ones(size(I)),v);

%% Compute levelset

hw=waitbar(0,'Computing levelset... ');
K=NaN(size(M));


logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);

sliceRange=find(~logicEmpty);
c=1;
numSteps=numel(sliceRange);
for qSlice=sliceRange    
    numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice    
    for qSub=1:1:numSubContours
        Vd=Vcs{qSlice}{qSub}; %Current contour
        if resampleOpt            
            d=pathLength(Vd);
            nResample=2*round(max(d(:))/min(v(:,[1 2])));
            [Vd]=evenlySampleCurve(V,nResample,'linear',1);
        end
        if ~isempty(Vd) %If it is not empty a levelset is computed 
            logicInside = inpolygon(x,y,Vd(:,1),Vd(:,2)); %Get logic for voxels inside contour           
            [D,~]=minDist([x(:) y(:)],Vd(:,[1 2])); %Compute distances to contour
            D(logicInside)=-D(logicInside); %Negate distances inside contour
            k_sub=reshape(D,[size(M,1) size(M,2)]); %The sub-levelset
            K(:,:,qSlice)=nanmin(K(:,:,qSlice),k_sub); %Add to current slice
        end
    end        
    waitbar(c/numSteps,hw,['Computing levelset. ',num2str(round((c/numSteps)*100)),'%']);    
    c=c+1;
end
close(hw);
