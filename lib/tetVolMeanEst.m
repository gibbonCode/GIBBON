function [regionA]=tetVolMeanEst(F,V)


% This function calculates the volume of an ideal regular tetrahedron with
% edge lengths (all equal) that match the mean edge lengths occuring for
% the input surface defined by F (faces) and V (vertices). 


%%
[edgeLengths]=patchEdgeLengths(F,V);
edgeLengthsMean=mean(edgeLengths);
meanProposedVolume=edgeLengthsMean^3./(6*sqrt(2)); %For a regular tetrahedron
regionA=meanProposedVolume;