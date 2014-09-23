function [M,S]=rigidbody_transform(X)

% function [M]=rigidbody_transform(V)
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
[U,S,V]=svd(X,0); %Singular value decomposition

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

R  = eye(4,4); R(1:3,1:3)=DCM;

%% Create translation rotation matrix
M  = T * R * eye(4,4);

end

