function F_cell=triLinearTri_F2D(VX,Vx,TRI)

F_cell=cell(size(TRI,1),1);
for qc=1:1:numel(F_cell);
    indNow=TRI(qc,:);
    X=VX(indNow,:); %Initial coordinates
    x=Vx(indNow,:); %Current coordinates
    F=triLinearTri_subF2D(X,x); %The deformation gradient tensor
    F_cell{qc}=F; %Store in cell array
end

end

function F=triLinearTri_subF2D(X,x)

% Define shape function set N=[1-r-s; r; s;];

%Hardcoded derivative
dN_dRST =[-1 -1;...
           1  0;...
           0  1];

% Compute derivatives of initial position vectors with respect to shape functions
dX_dRST=zeros(2,2);
for q=1:1:size(dN_dRST,1);
    dX_dRST=dX_dRST+(X(q,:)'*dN_dRST(q,:));
end

% Compute derivatives of shape functions with respect to initial position vectors
dN_dX=zeros(2,2);
for q=1:1:size(dN_dRST,1);
    dN_dX(q,:)=(dX_dRST'\dN_dRST(q,:)')';
end

%% DERIVE THE DEFORMATION GRADIENT TENSOR
F_2D=zeros(2,2);
for q=1:1:size(x,1)
    F_2D=F_2D+(x(q,:)'*dN_dX(q,:));
end
F=eye(3,3);
F(1:2,1:2)=F_2D;

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
