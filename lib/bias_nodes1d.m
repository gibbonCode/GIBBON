function [x]=bias_nodes1d(x,f)

% function [x]=bias_nodes1d(x,f)
% ------------------------------------------------------------------------
% This function biases the spacing for the entries in x using the factor f
% such that the output array goes from x(1) to x(end) but biased using x^f
% (i.e. the spacing decreasing depending on the magnitude of x). 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

min_x=min(x,[],2)*ones(1,size(x,2)); 
x=x-min_x; 
max_x=max(x,[],2)*ones(1,size(x,2)); 
x=max_x.*((x.^f)./(max_x.^f)); 
x=x+min_x;
 
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
