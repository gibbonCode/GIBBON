function [D]=logic2levelset(varargin)

% function [D]=logic2levelset(logicInside,voxelSize,voxelSizeResample)


%%

switch nargin
    case 1
        logicInside=varargin{1};
        voxelSize=[];
        voxelSizeResample=[];
    case 2
        logicInside=varargin{1};
        voxelSize=varargin{2};
        voxelSizeResample=[];
    case 3
        logicInside=varargin{1};
        voxelSize=varargin{2};
        voxelSizeResample=varargin{3};
end

if isempty(voxelSize)
    voxelSize=ones(1,3);
end

if isempty(voxelSizeResample)
    voxelSizeResample=mean(voxelSize);
end

if numel(voxelSizeResample)==1
    voxelSizeResample=voxelSizeResample*ones(1,3);
end

%%

siz=size(logicInside);

if ~max(voxelSize-mean(voxelSize))< max(eps(voxelSize))
    resampleOn=1;    
    [logicInside]=imageResample(logicInside,voxelSize,voxelSizeResample);
    scaleFactor=voxelSizeResample(1);
else
    resampleOn=0;
    scaleFactor=voxelSize(1);
end
logicInside=logicInside>0;

%Remove interior from logic
logicOn = bwmorph3(logicInside,'remove');

%Do distance transform on isotropic image
D = bwdist(logicOn,'euclidean');
D(logicInside)=-D(logicInside); %Negate distance inside
D=D.*scaleFactor; %Scale by voxel size

if resampleOn==1    
    sizNew=size(D);
    [J,I,K]=meshgrid(1:1:siz(2),1:1:siz(1),1:1:siz(3));
    [X,Y,Z]=im2cart(I,J,K,voxelSize);
    [I,J,K]=cart2im(X,Y,Z,voxelSizeResample);
    I=round(I);
    J=round(J);
    K=round(K);
    IND=reshape(sub2indn(sizNew,[I(:) J(:) K(:)],1),size(I));
    D=D(IND);
end

