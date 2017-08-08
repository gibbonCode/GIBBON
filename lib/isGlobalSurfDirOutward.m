function [L]=isGlobalSurfDirOutward(F,V)

% Not the best implementation at present. Vertices are offset allong the
% local normal direction by 1 10th of the smalles edge length. Then the
% volume before and after this operation is calculated. If the volume
% decreased the normals face the wrong way and the face orientation is thus
% flipped for smoothening. Contraction/inflation allong normal directions
% in this way does not always yield valid surfaces and hence volume
% computation may be inappropriate.

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
