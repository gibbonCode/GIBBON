function [V_norm]=vecnormalize(V)

if isvector(V)
    H=sqrt(sum(V.^2));
else
    H=sqrt(sum(V.^2,2));    
end
logicInvalid=H==0;

V_norm=V./H(:,ones(size(V,2),1));
V_norm(logicInvalid,:)=0; %Set invalid lengths to 0
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
