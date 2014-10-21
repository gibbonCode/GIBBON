function [v]=tetVolMeanEst(F,V)

% function [v]=tetVolMeanEst(F,V)
% ------------------------------------------------------------------------
%
% This function calculates the volume of an ideal regular tetrahedron with
% edge lengths (all equal) that match the mean edge lengths occuring for
% the input surface defined by F (faces) and V (vertices). 

%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/17
%------------------------------------------------------------------------
%%

[edgeLengths]=patchEdgeLengths(F,V);
edgeLengthsMean=mean(edgeLengths);
meanProposedVolume=edgeLengthsMean^3./(6*sqrt(2)); %For a regular tetrahedron
v=meanProposedVolume;