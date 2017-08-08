function cMap=grayJetc(varargin)

switch nargin
    case 0
        n=250;
    case 1
        n=varargin{1};
end

cMap=[1 0 0; 1 1 0; 0 1 0; 0 1 1; ];
cMap=[flipud(gray(size(cMap,1))); flipud(cMap)];
        

[cMap]=resampleColormap(cMap,n);
 
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
