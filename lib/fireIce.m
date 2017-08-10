function [cMap]=fireIce(varargin)

% function [cMap]=fireIce(n)
% ------------------------------------------------------------------------
% Creates the colormap data for n levels for the fire and ice colormap. Low
% values define a cold blue color while high values define a warm/hot
% color. The extrema of the colordata are the same color. Therefore this
% map is suitable for angular data where 360 may be equivalent to 0
% degrees.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/01/01
%------------------------------------------------------------------------

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

% cFire=[0 0 0; 1 0 0; 1 1 0; 1 1 1]; %Simple version
cFire=[0,0,0;0.112,0.0130,0;0.213,0.0240,0.00100;0.323,0.0370,0.00100;...
       0.421,0.0520,0.00100;0.526,0.0720,0.00200;0.620,0.100,0.00300;...
       0.693,0.131,0.0110;0.770,0.172,0.0190;0.837,0.217,0.0320;0.894,0.261,0.0430;...
       0.954,0.315,0.0520;0.993,0.379,0.0630;0.998,0.484,0.0780;1,0.574,0.110;...
       1,0.644,0.152;1,0.711,0.211;1,0.773,0.274;1,0.837,0.322;1,0.914,0.369;...
       1,0.953,0.443;1,0.982,0.539;0.999,1,0.646;0.985,1,0.772;1,1,1];
cMap=[rot90(cFire,2); cFire(2:end,:)];

[cMap]=resampleColormap(cMap,n);
 
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
