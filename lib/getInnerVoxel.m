function [indInternal]=getInnerVoxel(L,searchRadius,plotOn)

if any(L(:))
    %% Get an inner voxel
    
    %Construct search blurring kernel
    [x,y,z]=ndgrid(-searchRadius:1:searchRadius); %kernel coordinates
    d=sqrt(x.^2+y.^2+z.^2); %distance metric
    k=(d<=searchRadius); %Define kernel based on radius
    k=k./sum(k(:)); %Normalize kernel
    
    %Convolve interior image with kernel
    ML=convn(double(L),k,'same');
    ML(~L)=0;
    [~,indInternal]=max(ML(:)); %Kernel should yield max at "deep" (related to search radius) point
    
    %% Plotting image, mesh and inner point
    if plotOn==1
        voxelSize=[1 1 1];
        [I_in,J_in,K_in]=ind2sub(size(L),indInternal); %Convert to subscript coordinates
        [X_in,Y_in,Z_in]=im2cart(I_in,J_in,K_in,voxelSize);
        V_in=[X_in Y_in Z_in];
        
        %Plot settings        
        fontSize=20;
        markerSize1=50;
        faceAlpha1=0.75;
        
        cFigure        
        hold on;
        
        %Found voxel location
        hp=plotV(V_in,'r.','MarkerSize',markerSize1);
        
        %Local slices
        L_plot=false(size(L));
        L_plot(:,:,K_in)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=im2patch(L,L_plot,'sk');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        gpatch(Fm,Vm,Cm,'k',faceAlpha1);
        
        L_plot=false(size(L));
        L_plot(I_in,:,:)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=im2patch(L,L_plot,'si');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        gpatch(Fm,Vm,Cm,'k',faceAlpha1);
        
        L_plot=false(size(L));
        L_plot(:,J_in,:)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=im2patch(L,L_plot,'sj');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        gpatch(Fm,Vm,Cm,'k',faceAlpha1);
        
        colormap(gray(3));
        axisGeom(gca,fontSize);   
        camlight headlight;
        legend(hp,'Interior point','Location','NorthOutSide')
        drawnow;
    end
else
    indInternal=[];
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
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
