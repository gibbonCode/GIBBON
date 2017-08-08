function [I,J,K]=cart2im(X,Y,Z,v)

% function [I J K]=cart2im(X,Y,Z,v)
% ------------------------------------------------------------------------
% This function converts the cartesian coordinates X,Y,Z to image
% coordinates I,J,K using the voxel dimension v.
%
% X,Y,Z can be scalars, vectors or matrices. 
% v is a vector of length 3 where v(1), v(2) and v(3) correspond to the
% voxel dimensions in the x,y and z direction respectively. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2008/08/15
% ------------------------------------------------------------------------

I=(Y./v(1))+0.5;
J=(X./v(2))+0.5;
K=(Z./v(3))+0.5;
 
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
