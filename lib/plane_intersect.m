function X=plane_intersect(V1,V2,V3,N1,N2,N3)

%%

DET=ndet(N1,N2,N3);
DET(DET==0)=NaN;

%  X= (1./DET).*[(dot(V(:,1),N(:,1)).*cross(N(:,2),N(:,3))) + (dot(V(:,2),N(:,2)).*cross(N(:,3),N(:,1))) + (dot(V(:,3),N(:,3)).*cross(N(:,1),N(:,2)))]
 
X= ((1./DET)*ones(1,3)).* (...
     ((dot(V1,N1,2)*ones(1,3)).*cross(N2,N3,2)) + ...
     ((dot(V2,N2,2)*ones(1,3)).*cross(N3,N1,2)) + ...
     ((dot(V3,N3,2)*ones(1,3)).*cross(N1,N2,2))...
     );


end
 
 
 
%% <-- GIBBON footer text --> 
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
