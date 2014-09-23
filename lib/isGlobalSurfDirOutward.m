function [L]=isGlobalSurfDirOutward(F,V)

% Not the best implementation at present. Vertices are offset allong the
% local normal direction by 1 10th of the smalles edge length. Then the
% volume before and after this operation. If the volume decreased the
% normals face the wrong way and the face orientation is thus flipped for
% smoothening. Contraction/inflation allong normal directions in this way
% does not always yield valid surfaces and hence volume computation may be
% inappropriate.

%Compute edge lengths
[edgeLengths]=patchEdgeLengths(F,V);
minLength=min(edgeLengths(:));

%Get vertex normals
[~,~,N]=patchNormal(F,V);

%Compute volumes before and after "contraction/inflation" and flip faces if required.
[volFV1]=triSurfVolume(F,V); %Initial volume
[volFV2]=triSurfVolume(F,V+(minLength/10.*N)); %"contracted/inflated" volume

L=volFV2>volFV1; %if the volume increased the global direction is outward

end