function V=laplacian_smoothing(V,IND_V,L,n)

[I,J,v] = find(IND_V);

for i=1:n;
    Xp=accumarray({I,J},V(v,1),size(IND_V),[],NaN);
    Yp=accumarray({I,J},V(v,2),size(IND_V),[],NaN);
    Zp=accumarray({I,J},V(v,3),size(IND_V),[],NaN);
    Vp=[nanmean(Xp,2) nanmean(Yp,2) nanmean(Zp,2)];
    V=V+L.*(Vp-V);
end

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
