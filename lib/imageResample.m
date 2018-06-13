function [M_new]=imageResample(M,voxelSize,voxelSizeNew)

voxelSize=voxelSize(:)';
voxelSizeNew=voxelSizeNew(:)';

siz=size(M);
FOV=siz.*voxelSize; 
sizNew=ceil(FOV./voxelSizeNew);

[Jn,In,Kn]=meshgrid(1:1:sizNew(2),1:1:sizNew(1),1:1:sizNew(3));
[X,Y,Z]=im2cart(In(:),Jn(:),Kn(:),voxelSizeNew);
[I,J,K]=cart2im(X,Y,Z,voxelSize);
I=round(I);
J=round(J);
K=round(K);

logicValid=I<=siz(1) & J<=siz(2) & K<=siz(3);

I=I(logicValid);
J=J(logicValid);
K=K(logicValid);

IND=sub2ind(siz,I,J,K); 

M_new=nan(sizNew); 
M_new(logicValid)=M(IND);
