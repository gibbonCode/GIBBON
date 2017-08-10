function S=smoothHeaviside(varargin)

switch nargin
    case 1
        x=varargin{1};
        k=6;
        r=0;
    case 2
        x=varargin{1};
        k=varargin{2};
        r=0;
    case 3
        x=varargin{1};
        k=varargin{2};
        r=varargin{3};
end
k=k*2;
S=exp(2*k*(x-r))./(1+exp(2*k*(x-r)));

 
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
