function [V]=ellipseCoord(A,t)

% function [V]=ellipseCoord(A,t)
% ------------------------------------------------------------------------
% Calculates ellipse coordiantes for the angles in t based on the vector A
% which defines the centre coordinates, the radii and the angle
% respectively. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 2013/24/09
%------------------------------------------------------------------------

%%
x0=A(1);
y0=A(2);
x=A(3).*cos(t);
y=A(4).*sin(t);
V=[x(:) y(:) zeros(size(x(:)))];
[R,~]=euler2DCM([0 0 -A(5)]);
V=(R*V')';
V(:,1)=V(:,1)+x0;
V(:,2)=V(:,2)+y0;
 
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
