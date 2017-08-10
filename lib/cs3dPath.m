function [pp,t]=cs3dPath(V,p,w)

%Equivalent to CSCVN performance if p=1 and w=ones(size(V,1),1)

dt = sum((diff(V).^2).'); %Point spacing measure
t = cumsum([0,dt.^(1/4)]); %Curve length measure
% dt = sqrt(sum((diff(V).^2).')); %Point spacing measure
% t = cumsum([0,dt]); %Curve length measure
pp = csaps(t,V',p,[],w); %Smoothened ppform
 
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
