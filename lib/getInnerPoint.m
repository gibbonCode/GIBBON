function [varargout]=getInnerPoint(varargin)

% function [V_inner,M,G,ML]=getInnerPoint(F,V,searchRadius,voxelSize,plotOn)

%%

switch nargin
    case 2
        F=varargin{1};
        V=varargin{2};
        searchRadius=[];
        voxelSize=[];
        plotOn=0;
    case 3        
        F=varargin{1};
        V=varargin{2};
        searchRadius=varargin{3};
        voxelSize=[];
        plotOn=0;
    case 4            
        F=varargin{1};
        V=varargin{2};
        searchRadius=varargin{3};
        voxelSize=varargin{4};
        plotOn=0;
    case 5
        F=varargin{1};
        V=varargin{2};
        searchRadius=varargin{3};
        voxelSize=varargin{4};
        plotOn=varargin{5};
    otherwise
        error('Wrong number of input arguments');
end

if isa(F,'cell')
    [F,V,C]=joinElementSets(F,V);
else
    C=ones(size(F,1),1);
end

if isempty(voxelSize)
   D=patchEdgeLengths(F,V); 
   voxelSize=mean(D)/2;
end

if isempty(searchRadius)
   searchRadius=voxelSize*3;
end

%% Get an inner point

%Convert patch data to image
[M,G,~]=patch2Im(F,V,C,voxelSize);
L_test=(M==1); %Logic for interior voxels of surface mesh
imOrigin=G.origin; %Image origin can be used to allign image with surface

%Construct search blurring kernel
[x,y,z]=ndgrid(-searchRadius:voxelSize:searchRadius); %kernel coordinates
d=sqrt(x.^2+y.^2+z.^2); %distance metric
k=(d<=searchRadius); %Define kernel based on radius
k=k./sum(k(:)); %Normalize kernel

%Convolve interior image with kernel
ML=convn(double(L_test),k,'same');

% Alternative based on boundary distance transform (slower)
%ML = bwdist(M==0,'euclidean'); ML(~L_test)=0;

ML(M~=1)=NaN; %Set other sites to NaN
[~,indInternal]=gnanmax(ML(:)); %Kernel should yield max at "deep" (related to search radius) point
[I_in,J_in,K_in]=ind2sub(size(M),indInternal); %Convert to subscript coordinates

%Derive point coordinates
V_inner=zeros(1,3);
[V_inner(1),V_inner(2),V_inner(3)]=im2cart(I_in,J_in,K_in,voxelSize*ones(1,3));
V_inner=V_inner+imOrigin(ones(size(V_inner,1),1),:);

%% Plotting image, mesh and inner point
if plotOn==1    
    
    fontSize=20;
    markerSize1=50;%round(voxelSize*25);
    faceAlpha1=0.75;
    faceAlpha2=0.5;
    
    cFigure;    
    hold on;
    
    gpatch(F,V,0.5*ones(1,3),'none',faceAlpha2);
    plotV(V_inner,'k.','MarkerSize',markerSize1);
    
    L_plot=false(size(ML));
    L_plot(:,:,K_in)=1;
    L_plot=L_plot & ~isnan(ML);
    [Fm,Vm,Cm]=ind2patch(L_plot,double(ML),'sk');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    
    gpatch(Fm,Vm,Cm,'k',faceAlpha1);    
    
    L_plot=false(size(ML));
    L_plot(I_in,:,:)=1;
    L_plot=L_plot & ~isnan(ML);
    [Fm,Vm,Cm]=ind2patch(L_plot,ML,'si');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    gpatch(Fm,Vm,Cm,'k',faceAlpha1);    
    
    L_plot=false(size(ML));
    L_plot(:,J_in,:)=1;
    L_plot=L_plot & ~isnan(ML);
    [Fm,Vm,Cm]=ind2patch(L_plot,ML,'sj');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    gpatch(Fm,Vm,Cm,'k',faceAlpha1);    
    
    colormap(gjet(250)); colorbar;     
    axisGeom(gca,fontSize);
    drawnow;    
end

varargout{1}=V_inner;
varargout{2}=M;
varargout{3}=G;
varargout{4}=ML;
varargout{5}=voxelSize;
 
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
