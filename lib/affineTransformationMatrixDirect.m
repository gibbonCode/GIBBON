function [T]=affineTransformationMatrixDirect(V1,V2)

%%
%Force input to 3D
if size(V1,2)==2
    V1(:,3)=0; 
end

if size(V2,2)==2
    V2(:,3)=0; 
end

%Expand to nx4
V1_M=V1;
V1_M(:,4)=1; 
V2_M=V2;
V2_M(:,4)=1; 

%Get transformation using left devide
T=(V1_M\V2_M)';

 
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
