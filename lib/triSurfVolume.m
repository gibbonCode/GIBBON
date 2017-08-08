function [VV]=triSurfVolume(F,V)

% function [VV]=triSurfVolume(F,V)
% ------------------------------------------------------------------------
%Volume derivation based on Gauss divergence theorem
%
%%% EXAMPLE
%
% r=3; %sphere radius
% n=2; %Refinements   
% [F,V,~]=geoSphere(4,r);
% 
% VVt=(4/3)*pi*r.^3; %Theoretical volume of the sphere
% [VV]=triSurfVolume(F,V); %estimate based on triangulated surface
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 26/11/2013
%------------------------------------------------------------------------
%%

[N,~]=trinorm(F,V); %Face normals
aa=tri_area(F,V); %Areas
Zm = (V(F(:,1),3)+V(F(:,2),3)+V(F(:,3),3))./3; %Mean Z
Nz = N(:,3); %Z component of normal
vv = aa.*Zm.*Nz; %Contributions
VV = sum(abs(vv)); %Total volume

 
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
