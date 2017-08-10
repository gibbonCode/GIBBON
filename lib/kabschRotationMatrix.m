function [Q]=kabschRotationMatrix(V1,V2)

%2017/01/19 Fixed bug in relation to forcing right handed coordinate
%system. This will avoid inverting as well. 

%%
%Force input to 3D
if size(V1,2)==2
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

%% 

%Force input to 3D
if size(V1,2)==2
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

%Centering on the mean
V1_m=mean(V1,1);
V1=V1-V1_m(ones(size(V1,1),1),:);

V2_m=mean(V2,1);
V2=V2-V2_m(ones(size(V2,1),1),:);

%% Use Kabsch algorithm

A=V1'*V2;

[U,~,V] = svd(A);

d=sign(det(U'*V));
D=eye(3,3);
D(3,3)=d; 
Q=V*D*U';

 
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
