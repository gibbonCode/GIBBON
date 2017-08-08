function [F_fixed,L_fixed]=fixNormalsOutward(F,V,fixOpt)

%Assumes F and V describe a shape that can be appropirately defined using a
%spherical or polar coordinate system, as indicated by the choice fixOpt.
%For instance a closed spherical shape for open cylindrical shape. 
%
% The function uses the fact that the radius of the normal vector origin
% should be smaller than the radius of the normal vector tip to determine
% whether or not a normal vector needs to be flipped. 

%Derive current face normals
[N,V_starts]=patchNormal(F,V);

%Vector lengths may need to be scaled, scaling factor should be
%arbitrary if and 100 should work is shape is sufficiently
%spherical. For severly "deformed" or noisy shapes this may
%have to be smaller, but one could argue that then the
%spherical coordinate system mapping is not appropriate.
f=100;
switch fixOpt
    case 'p' %polar
        R_starts = sqrt(V_starts(:,1).^2 + V_starts(:,2).^2); %Start radii
        N=(min(R_starts(:))/f)*N; %Scaling normal lengths
        V_ends=V_starts+N; %vector end points
        R_ends = sqrt(V_ends(:,1).^2 + V_ends(:,2).^2); %End radii        
    case 's' %spherical         
        R_starts = sqrt(V_starts(:,1).^2 + V_starts(:,2).^2 + V_starts(:,3).^2); %Start radii
        N=(min(R_starts(:))/f)*N; %Scaling normal lengths
        V_ends=V_starts+N; %vector end points
        R_ends = sqrt(V_ends(:,1).^2 + V_ends(:,2).^2 + V_ends(:,3).^2); %End radii
end
L_fixed=R_ends<R_starts;
F_fixed=F;
F_fixed(L_fixed,:)=fliplr(F(L_fixed,:)); %Invert faces whose normal vector end radii are smaller than their origin

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
