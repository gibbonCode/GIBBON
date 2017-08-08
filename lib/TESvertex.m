function [TES_dist]=TESvertex(TES,x,y,z)

x_vert=x(TES);
y_vert=y(TES);
z_vert=z(TES);

c_from=1:1:size(TES,2);
c_upto=[c_from(end) c_from(1:end-1)];
TES_dist=zeros(size(TES,1),numel(c_from));
for c=c_from;
    dc=sqrt((x_vert(:,c_from(c))-x_vert(:,c_upto(c))).^2+(y_vert(:,c_from(c))-y_vert(:,c_upto(c))).^2+(z_vert(:,c_from(c))-z_vert(:,c_upto(c))).^2);
    TES_dist(:,c)=dc;
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
