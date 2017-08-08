function [L]=isPolyClockwise(V)

mean_V=mean(V,1);
V=V-mean_V(ones(size(V,1),1),:);

%Compute angles
T = atan2(V(:,2),V(:,1));

%Unwrap
T = unwrap(T);

%Is the mean of the derivative smaller than 0?
L=mean(diff(T))<0;
 
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
