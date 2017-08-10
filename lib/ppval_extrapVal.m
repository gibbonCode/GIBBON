function v=ppval_extrapVal(pp,xx,interpLim,extrapVal)

%Evaluate piecewise polynomial
v=ppval(pp,xx);

%Set extrapolation values 
v(xx<interpLim(1))=extrapVal(1);
v(xx>interpLim(2))=extrapVal(2);
 
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
