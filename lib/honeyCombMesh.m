function [Fh,Vh]=honeyCombMesh(minV,maxV,pointSpacing)

%% CREATE TRIANGULATION

maxV(2)=maxV(2)+pointSpacing; 
[F,V]=triMeshEquilateral(minV,maxV,pointSpacing);

%% GET DUAL FOR HONEY-COMB

[Vh,Fd]=patch_dual(V,F);
numVert=cellfun(@(x) size(x,2),Fd);
Fh=Fd{numVert==6};

