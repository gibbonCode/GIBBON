function [Vm]=rigidTransformVertices(V,T,v)

Vm=V; 
[I,J,K]=cart2im(V(:,1),V(:,2),V(:,3),v); %Convert to image coordinates
IJK=[I(:) J(:) K(:) ones(size(I(:)))]; %Prepare for mapping
IJK_mapped=(T*IJK')'; %Do mapping
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(IJK_mapped(:,1),IJK_mapped(:,2),IJK_mapped(:,3),v); %Convert mapped image coordinates back to "Cartesian" coordinates
 
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
