function axisGeom(varargin)

switch nargin
    case 0
        h=gca;
        fontSize=15;
    case 1
        h=varargin{1};
        fontSize=15;
    case 2
        h=varargin{1};
        fontSize=varargin{2};
end

if isempty(h)
    h=gca;
end

axes(h);
view(3);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
set(gca,'FontSize',fontSize);
axis equal; axis vis3d; axis tight;
grid on; box on; 
 
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
