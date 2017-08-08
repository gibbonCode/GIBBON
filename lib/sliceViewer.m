function hf=sliceViewer(varargin)


%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
        viewerType=1;
    case 2
        M=varargin{1};
        v=varargin{2};
        viewerType=1;
    case 3
        M=varargin{1};
        v=varargin{2};
        viewerType=varargin{3};
end
        
%% Start viewer

switch viewerType
    case 1
        hf=sv2(M,v);
    case 2
        hf=sv3(M,v);
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
