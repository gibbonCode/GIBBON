function [N,Vn]=trinorm(F,V)

%N.B. if F does not describe triangles this functions uses first three
%vertices as a triangular description


%Getting triangle surface normal (cross product of two edge vectors)
vec1=[V(F(:,2),1)-V(F(:,1),1)  V(F(:,2),2)-V(F(:,1),2)  V(F(:,2),3)-V(F(:,1),3)];
vec2=[V(F(:,3),1)-V(F(:,1),1)  V(F(:,3),2)-V(F(:,1),2)  V(F(:,3),3)-V(F(:,1),3)];
N=cross(vec1,vec2,2);

%Normalizing vector length
N=N./(sqrt(sum(N.^2,2))*ones(1,size(N,2)));

%Midface coordinates for normal vectors (mean of each face)
X=V(:,1); Y=V(:,2); Z=V(:,3);
Vn=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];

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
