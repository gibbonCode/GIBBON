function scf(h)

%function scf(h)
%SCF Get handle to current figure.
%   SCF(h) Makes the figure specified by the handle h the current figure
%   (without making it appear). 
%
%   See also FIGURE, CLOSE, GCF, CLF, GCA, GCBO, GCO, GCBF.


set(0,'CurrentFigure',h);
 
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
