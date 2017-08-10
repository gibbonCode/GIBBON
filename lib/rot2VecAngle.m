function [theta,w]=rot2VecAngle(R)

theta=acos(0.5*(trace(R)-1));
theta=real(theta);
w=(1./(2*sin(theta)))*([R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]);
w=vecnormalize(w);
 
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
