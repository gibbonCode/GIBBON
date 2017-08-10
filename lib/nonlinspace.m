function y=nonlinspace(range_in,range_out,fform,n)

x=linspace(range_in(1),range_in(2),n); 
eval(['y=',fform,';']);

y=y-y(1);
y=range_out(1)+(abs(range_out(1)-range_out(2))*(y./y(end)));

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
