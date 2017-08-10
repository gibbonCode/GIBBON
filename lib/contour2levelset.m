function [K]=contour2levelset(M,v,Vcs,levelSetType)


%%
resampleOpt=0;
siz=size(M); %Size of the image dataset

switch levelSetType
    case 1 %Slice-by-slice 2D distances
        %% Get image coordinates
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
                if ~isempty(Vd) %If it is not empty a levelset is computed
                    if resampleOpt
                        d=pathLength(Vd);
                        nResample=2*round(max(d(:))/min(v(:,[1 2])));
                        [Vd]=evenlySampleCurve(V,nResample,'linear',1);
                    end                    
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
    case 2 %3D distances       
%         %% Get image coordinates
%         [J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));
%         
%         %Convert to cartesian coordinates using voxel size if provided
%         [x,y,z]=im2cart(I,J,K,v);
%         
        %% Compute levelset
        
        [logicInside]=contour2logic(M,v,Vcs);
        
        logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);        
        sliceRange=find(~logicEmpty);
        
        VD=[];
        for qSlice=sliceRange
            numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice
            for qSub=1:1:numSubContours
                Vd=Vcs{qSlice}{qSub}; %Current contour
                if ~isempty(Vd) %If it is not empty a levelset is computed
                    if resampleOpt
                        d=pathLength(Vd);
                        nResample=2*round(max(d(:))/min(v(:,[1 2])));
                        [Vd]=evenlySampleCurve(V,nResample,'linear',1);
                    end                                        
                    VD=[VD;Vd];                    
                end
            end        
        end
        
        %Use distance transform
        D1 = bwdist(logicInside,'euclidean');                
        D2 = bwdist(~logicInside,'euclidean');
        K=zeros(siz);
        K(logicInside)=-D2(logicInside);
        K(~logicInside)=D1(~logicInside);
        
%         [D,~]=minDist([x(:) y(:) z(:)],VD); %Compute distances to contour
%         D(logicInside)=-D(logicInside); %Negate distances inside contour
%         K=reshape(D,siz); %The sub-levelset        
        
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
