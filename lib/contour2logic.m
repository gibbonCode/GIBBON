function [L]=contour2logic(M,v,Vcs)


%%

siz=size(M); %Size of the image dataset

%Get image coordinates
[J,I]=meshgrid(1:1:siz(2),1:1:siz(1));

%Convert to cartesian coordinates using voxel size if provided
[x,y,~]=im2cart(I,J,ones(size(I)),v);

%% Compute levelset

hw=waitbar(0,'Computing logic... ');
L=false(size(M));

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
            L(:,:,qSlice)=L(:,:,qSlice) | logicInside; %Add to current slice
        end
    end        
    waitbar(c/numSteps,hw,['Computing logic. ',num2str(round((c/numSteps)*100)),'%']);    
    c=c+1;
end
close(hw);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
