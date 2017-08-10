function [D]=distND(V1,V2)

D=zeros(size(V1,1),size(V2,1));

for q=1:1:size(V1,2) %For all dimensions
    A=V1(:,q);
    B=V2(:,q);
    
    AB=A(:,ones(1,size(B,1)));
    BA=B(:,ones(1,size(A,1)))';
    
    D=D+(AB-BA).^2;
end
D=sqrt(D);
 
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
