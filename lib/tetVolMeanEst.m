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
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
