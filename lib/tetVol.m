function VE=tetVol(E,V)

% function VE=tetVol(E,V)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the tetrahedral elements specified by the
% element matrix E and the vertices V.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
%------------------------------------------------------------------------

X=V(:,1); Y=V(:,2); Z=V(:,3);
XE=X(E); YE=Y(E); ZE=Z(E);
if size(E,1)==1 %Transpose in this special case
   XE=XE'; YE=YE'; ZE=ZE';
end
A=[XE(:,1) YE(:,1) ZE(:,1)];
B=[XE(:,2) YE(:,2) ZE(:,2)];
C=[XE(:,3) YE(:,3) ZE(:,3)];
D=[XE(:,4) YE(:,4) ZE(:,4)];

VE=abs(dot((A-D),cross((B-D),(C-D),2),2))./6;
 
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
