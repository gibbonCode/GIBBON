function [K]=contour2levelset(M,voxelSize,Vcs,levelSetType)

%-------------------------------------------------------------------------
% function [K]=contour2levelset(M,v,Vcs,levelSetType)
% This function converts contours to a levelset image. The levelset image
% represents distances to the contour. 
%
% The input consists of: 
% * the 3D image M
% * The vozelSize a 1x3 vector specifying the size of the voxels in the
% row, column, and slice direction.
% * A cell array Vcs containing one or more contours per slice. If the
% image has n slices then Vcs should be an nx1 cell array, i.e. contours
% are defined for each slice. 
% * The levelSetType, this variable determines the levelset type to use. If
% it equals 1 then the levelset image is formulated using 2D distances for
% each slice. If it equals 2 3D distances are computed (most
% computationally intensive). If it equals 3 the distance metric is the
% minimum between a 3D distance is based on distance transform and the 2D
% per slice distances. 
%
% Change log: 
% 2019/08/26 Added description to function. Updated handling of input M,
% can now instead specify the size of M too. 
%-------------------------------------------------------------------------


%%
resampleOpt=1;

%Get size of the image
if isvector(M) %Assume M is really size(M)
    siz=M;
else %Compute size(M) from M
    siz=size(M); 
end

switch levelSetType
    case 1 %Slice-by-slice 2D distances
        
        % Get image coordinates
        [J,I]=meshgrid(1:1:siz(2),1:1:siz(1));
        
        %Convert to cartesian coordinates using voxel size if provided
        [x,y,~]=im2cart(I,J,ones(size(I)),voxelSize);
        
        %Compute levelset
        hw=waitbar(0,'Computing distances... ');
        K=NaN(siz);
        
        logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);
        
        sliceRange=find(~logicEmpty);
        c=1;
        numSteps=numel(sliceRange);
        for qSlice=sliceRange
            numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice
            for qSub=1:1:numSubContours
                Vd=Vcs{qSlice}{qSub}; %Current contour
                if ~isempty(Vd) %If it is not empty a levelset is computed
                    logicInside = inpolygon(x,y,Vd(:,1),Vd(:,2)); %Get logic for voxels inside contour
                    if resampleOpt==1
                        d=pathLength(Vd);
                        nResample=2*round(max(d(:))/min(voxelSize([1 2])));
                        [Vd]=evenlySampleCurve(Vd,nResample,'linear',1);
                    end
                    [D,~]=minDist([x(:) y(:)],Vd(:,[1 2])); %Compute distances to contour
                    D(logicInside)=-D(logicInside); %Negate distances inside contour
                    k_sub=reshape(D,[size(M,1) size(M,2)]); %The sub-levelset
                    K(:,:,qSlice)=gnanmin(K(:,:,qSlice),k_sub); %Add to current slice
                end
            end
            waitbar(c/numSteps,hw,['Computing distances. ',num2str(round((c/numSteps)*100)),'%']);
            c=c+1;
        end
        close(hw);
    case 2 %3D Euclidean distance transform on logic for voxels inside contours
        [logicIn,logicOn]=contour2logic(M,voxelSize,Vcs); %Compute logic to denote inside
        logicInside=logicIn | logicOn;
        [K]=logic2levelset(logicInside,voxelSize,min(voxelSize)); %Compute distance transform on logic
    case 3 %A hybrid between method 1 and the distance transform
        %Compute logic to denote inside
        [logicIn,logicOn]=contour2logic(M,voxelSize,Vcs);
        logicInside=logicIn | logicOn;
        
        %Compute distance transform on logic as first estimate
        K=logic2levelset(logicInside,voxelSize,min(voxelSize));
        K=abs(K); %Force absolute so minimum is easy
        
        %Get boundary voxels
        try
            logicBoundary = bwmorph3(logicInside,'remove'); %New in R2018a
        catch
            logicBoundary=logicRemoveInterior(logicInside); %GIBBON alterative
        end
        
        %Grow boundary logic once using convolution
        hg=zeros(3,3,3); hg([5 11 13 14 15 17 23])=1; %3D mask
        logicBoundary=convn(double(logicBoundary),hg,'same')>0; %grow
        
        % Get image coordinates
        [JJ,II,KK]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));
        
        %Convert to cartesian coordinates using voxel size if provided
        [x,y,z]=im2cart(II,JJ,KK,voxelSize);
        
        %Collect contour points
        logicEmpty=cellfun(@(x) isempty(x{1}),Vcs);
        sliceRange=find(~logicEmpty);
        D=nan(size(K));
        hw=waitbar(0,'Computing distances... ');
        c=1; numSteps=numel(sliceRange);
        for qSlice=sliceRange
            numSubContours=numel(Vcs{qSlice}); %Number of sub-contours for the current slice
            for qSub=1:1:numSubContours
                Vd=Vcs{qSlice}{qSub}; %Current contour
                if ~isempty(Vd) %If it is not empty a levelset is computed
                    if resampleOpt==1
                        d=pathLength(Vd);
                        nResample=2*round(max(d(:))/min(voxelSize([1 2])));
                        [Vd]=evenlySampleCurve(Vd,nResample,'linear',1);
                    end
                    if qSlice==sliceRange(1) || qSlice==sliceRange(end)
                        [~,Vdd]=regionTriMesh2D({Vd(:,[1 2])},min(voxelSize([1 2])),1,0);
                        Vdd(:,3)=mean(Vd(:,3));
                        Vd=Vdd;
                    end
                    [D_set,~]=minDist([x(logicBoundary) y(logicBoundary) z(logicBoundary)],Vd); %Compute distances to contour
                    D(logicBoundary)=min(D(logicBoundary),D_set);
                end
            end
            waitbar(c/numSteps,hw,['Computing distances. ',num2str(round((c/numSteps)*100)),'%']);
            c=c+1;
        end
        close(hw);
        K(logicBoundary)=D(logicBoundary);
        K(logicInside)=-K(logicInside); %Negate distances inside contour
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
