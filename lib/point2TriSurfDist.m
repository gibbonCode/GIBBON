function D=point2TriSurfDist(F1,V1,V2)

%Get face normals
[N1,~]=trinorm(F1,V1);

%Compute face centres
X=V1(:,1); Y=V1(:,2); Z=V1(:,3);
XF=X(F1); YF=Y(F1); ZF=Z(F1);

%Position vectors for face centres
v1=[mean(XF,2) mean(YF,2) mean(ZF,2)];

%Find closest triangles
Dm=dist(V2,v1');
[~,indMin]=min(Dm,[],2); %Get indices of closest

%Difference vectors between face centres and points
x2=V2-v1(indMin,:);

%Find distance to closest triangles
D=abs(dot(N1(indMin,:),x2,2));
 
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
