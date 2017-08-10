function [M,S]=rigidbody_transform(X)

% function [M]=rigidbody_transform(X)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 22/04/2011
% ------------------------------------------------------------------------

%% Determine translation matrix
OR=mean(X,1);

%Defining translation matrix
T  = [1 0 0 OR(1);...
      0 1 0 OR(2);...
      0 0 1 OR(3);...
      0 0 0 1];

%% Determine rotation matrix

X=[X(:,1)-OR(1) X(:,2)-OR(2) X(:,3)-OR(3)]; %Centre points around mean
[~,S,V]=svd(X,0); %Singular value decomposition

%Defining direction cosine matrix

if 1-V(3,3)<eps('double')
    DCM=eye(3,3);
else
    rz=V(:,3); rz=rz./sqrt(sum(rz.^2)); %surface normal
    r=V(:,2); r=r./sqrt(sum(r.^2));
    rx=cross(rz,r);rx=rx./sqrt(sum(rx.^2));
    ry=cross(rx,rz);ry=ry./sqrt(sum(ry.^2));
    DCM=[rx(:) ry(:) rz(:)];
end

% N=-V(:,3)./V(3,3); %Surface normal
% rx=[1 0 N(1)]; rx=rx./sqrt(sum(rx.^2));
% ry=[0 1 N(2)]; ry=ry./sqrt(sum(ry.^2));
% rz=cross(rx,ry); rz=rz./sqrt(sum(rz.^2));
% DCM=[rx(:) ry(:) rz(:)];

R  = eye(4,4); 
R(1:3,1:3)=DCM;

%% Create translation rotation matrix
M  = T * R * eye(4,4);

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
