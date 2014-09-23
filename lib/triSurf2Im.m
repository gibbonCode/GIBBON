function [varargout]=triSurf2Im(varargin)

% function [M,G,bwLabels]=triSurf2Im(F,V,voxelSize,imOrigin,imSiz)
% -----------------------------------------------------------------------
% This function converts the input triangulated surface, specified by the
% faces F and the vertices V into an image. The image size, origin and
% voxel size can be specified. If not specified the image voxel size is set
% to the mean edge length, the size is set to fit the object (with an extra
% voxel spacing all around) and the origin is choosen in the corner of the
% tightly fit image.  
% If required the surface is resampled to subvoxel resolution using the
% |subtri| function. Then surface vertices are simply mapped to an image
% coordinate system and since they are densely sampled with respect to the
% voxel size form an enclosing boundary of voxels. The exterior, boundary
% and interior are then formulated using the bwlabeln function. The
% exterior, boundary and intertior voxels are labelled using 0's, 1's and
% 2's in the output image M. The second output G is a structure containing
% the fields G.voxelSize and G.origin which form the geometric information
% for the image. 
% 
% 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 28/08/2013
%------------------------------------------------------------------------

if nargin<2 || nargin>5
        error('Wrong number of input arguments!');
end

%Get faces and vertices
F=varargin{1};
V=varargin{2};

%Removing unused points
[F,V]=removeNotIndexed(F,V);

%Checking edge lenghts of surface
[edgeLengths]=patchEdgeLengths(F,V);

maxEdgeLength=max(edgeLengths(:));
meanEdgeLength=mean(edgeLengths(:));

if nargin >2
    voxelSize=varargin{3};
else
    voxelSize=meanEdgeLength;
end

if nargin >3
    imOrigin=varargin{4};
else
    %Determine surface set coordinate minima
    minV=min(V,[],1);
    
    %Determine shift so all coordinates are positive
    imOrigin=(minV-voxelSize);
end

    
%%

%Resample surface if voxelsize is small with respect to edgelenghts
n=maxEdgeLength/voxelSize;
if (n-1)>eps(1)
    n=ceil(n);
    [~,V]=subtri(F,V,n);
end

%Shift points using origin
V=V-imOrigin(ones(size(V,1),1),:);

%Convert to image coordinates
V_IJK=V;
[V_IJK(:,1),V_IJK(:,2),V_IJK(:,3)]=cart2im(V(:,1),V(:,2),V(:,3),voxelSize*ones(1,3));

%Rounding image coordinates to snap to voxel
V_IJK=round(V_IJK);

%Determine image size if not provided
if nargin<5    
    siz=max(V_IJK,[],1)+1;
else
    siz=varargin{5};
    
    %Remove invalid indices
    V_IJK=V_IJK(V_IJK(:,1)<=siz(1) & V_IJK(:,1)>=1,:);
    V_IJK=V_IJK(V_IJK(:,2)<=siz(2) & V_IJK(:,2)>=1,:);
    V_IJK=V_IJK(V_IJK(:,3)<=siz(3) & V_IJK(:,3)>=1,:);       
end

%Get linear indices of points
indVertices=sub2ind(siz,V_IJK(:,1),V_IJK(:,2),V_IJK(:,3));

%Create surface boundary image
logicVertices=false(siz);
logicVertices(indVertices)=1;

%Create boundary, interior and exterior image  
labeledImage = bwlabeln(~logicVertices,6); %Get labels for ~L image which will segment interior, exterior and boundary
uniqueLabels=unique(labeledImage(:));

labelsBoundary=labeledImage(logicVertices); %The label numbers for the boundary

indExteriorVoxel=1; %First is outside since image is at least a voxel too big on all sides
labelExterior=labeledImage(indExteriorVoxel); %The label number for the exterior

labelsInterior=uniqueLabels(~ismember(uniqueLabels,[labelsBoundary(:); labelExterior])); %Labels for the interior (possibly multiple)

M=zeros(size(logicVertices)); %The exterior is set to 0
M(logicVertices)=1; %The boundary is set to 1
M(ismember(labeledImage,labelsInterior))=2; %Interior is set to 2

%Storing image geometry metrics
G.voxelSize=voxelSize;
G.origin=imOrigin;

%Create output
labeledImage(M==0)=NaN;
varargout{1}=M;
varargout{2}=G;
varargout{3}=labeledImage;
    