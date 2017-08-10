function P=inv_parbound(Pb,Ps)

w=abs(Ps.ub-Ps.lb); %Interval width
b=Ps.f*w; %Point at which the output interval width equals t*w
a=tan(Ps.t*0.5*pi);
P=Ps.c+((b./a)*tan((pi*(Pb-Ps.lb)./w)-(pi/2)));




 
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
