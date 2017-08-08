function [Z]=rozenbrock(X,Y)

% -----------------------------------------------------------------------
% function [Z]=rozenbrock(X,Y)
%
% "In mathematical optimization, the Rosenbrock function is a non-convex
% function used as a test problem for optimization algorithms. It is also
% known as Rosenbrock's valley or Rosenbrock's banana function. This
% function is often used to test performance of optimization algorithms.
% The global minimum is inside a long, narrow, parabolic shaped flat
% valley. To find the valley is trivial, however to converge to the global
% minimum is difficult.It is defined by Z=(1-X.^2)+100.*((Y-(X.^2)).^2)." 
% 
% From: http://en.wikipedia.org/wiki/Rosenbrock_function
%
% Kevin Moerman
% kevinmoerman@hotmail.com
% -----------------------------------------------------------------------

Z=(1-X.^2)+100.*((Y-(X.^2)).^2);
 
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
