function [Dt]=patchVectorTangent(F,V,D,N)

%Derive tangent contribution through dot-product with normal vectors

if isempty(N)
    [~,~,N]=patchNormal(F,V); %Get current vertex normals    
end

Dn_mag=dot(D,N,2); %Allong normal displacement magnitudes
Dn=Dn_mag(:,ones(1,3)).*N; %Normal direction displacement vectors
Dt=D-Dn; %Tranverse or tangential only displacement vectors
 
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
