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

 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
