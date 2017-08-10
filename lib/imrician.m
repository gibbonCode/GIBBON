function M_rice=imrician(M,s)

% function M_rice=imrician(M,s)
% ------------------------------------------------------------------------
% IMCRICIAN Random samples from the Rice/Rician probability distribution.
%
% R ~ Rice(v, s) if R = sqrt(X^2 + Y^2), where X ~ N(v*cos(a), s^2) and
% Y ~ N(v*sin(a), s^2) are independent normal distributions (any real a).
%
% Reference: http://en.wikipedia.org/wiki/Rice_distribution
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 15/04/2009
% ------------------------------------------------------------------------

X_GAUSS=s.*randn(size(M))+M;
Y_GAUSS=s.*randn(size(M));
M_rice = hypot(X_GAUSS,Y_GAUSS); 
 
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
