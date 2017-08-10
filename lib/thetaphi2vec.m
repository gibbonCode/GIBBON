function [r,c]=thetaphi2vec(THETA,PHI)

%--------------------------------------------------------------------------
% function [vx,vy]=thetaphi2vec(THETA,PHI)
% 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/02/2010
%-------------------------------------------------------------------------- 


%%

r = [cos(THETA).*cos(PHI)  -sin(THETA) cos(THETA).*sin(PHI)];
c = [sin(THETA).*cos(PHI)  cos(THETA)  sin(THETA).*sin(PHI)];

end
 
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
