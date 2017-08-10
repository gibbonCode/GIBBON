function [varargout]=rigidTransformationMatrixDirect(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2 
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

[Q]=kabschRotationMatrix(V1,V2);

V1_m=mean(V1,1);
V2_m=mean(V2,1);

T1=eye(4,4);
T1(1:3,end)=V2_m(:);

R=eye(4,4);
R(1:3,1:3)=Q;

T2=eye(4,4);
T2(1:3,end)=-V1_m;

T=T1*R*T2;

varargout{1}=T;
varargout{2}=R(1:3,1:3);

 
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
