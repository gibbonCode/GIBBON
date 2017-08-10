function A=tri_area(F,V)

if size(V,2)==2
    V(:,3)=0;
end

V12=V(F(:,2),:)-V(F(:,1),:);
V13=V(F(:,3),:)-V(F(:,1),:);
A=0.5.*sqrt(sum(cross(V12,V13,2).^2,2));
 
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
