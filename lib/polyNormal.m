function [varargout]=polyNormal(V_poly)

% [N,Vn,Nv]=polyNormal(V_poly)
% ------------------------------------------------------------------------
% Normals directions are derived based on cross product of curve segments
% with the outware z-direction. A horizontal and clockwise curve segment
% therefore results in an upward normal. The normals are provided for each
% segment or for each point on the curve. 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2017/03/16 Created
%------------------------------------------------------------------------

%%

zDir=[0 0 1];

if size(V_poly,2)==2 %Cope with 2D
    V_poly(:,3)=0;
end

[Vn,U]=polyDerivative(V_poly,1);
N=cross(zDir(ones(size(U,1),1),:),U);
N=vecnormalize(N);

varargout{1}=N(:,[1 2]); 
varargout{2}=Vn(:,[1 2]); 

if nargout==3
    [~,U]=polyDerivative(V_poly,3);
    Nv=cross(zDir(ones(size(U,1),1),:),U);
    Nv=vecnormalize(Nv);
    varargout{3}=Nv(:,[1 2]);
end

end

function [V,U]=polyDerivative(Vg,dirOpt)

switch dirOpt
    case 1 %Forward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=diff(Vg,1,1);        
    case 2 %Backward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=flipud(diff(flipud(Vg),1,1));        
    case 3 %Central
        V=Vg;
        Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
        Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
        U=Uf;
        U(:,:,2)=Ub;
        U=nanmean(U,3);        
end
U=vecnormalize(U);

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
