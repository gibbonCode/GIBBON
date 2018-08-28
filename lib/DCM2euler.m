function [a]=DCM2euler(R)


%%

% cosTheta = cos(a);
% sinTheta = sin(a);
%         
% R(1,1,:) = cosTheta(:,2).*cosTheta(:,3);
% R(1,2,:) = -cosTheta(:,2).*sinTheta(:,3);
% R(1,3,:) = sinTheta(:,2);
% R(2,1,:) = cosTheta(:,1).*sinTheta(:,3) + cosTheta(:,3).*sinTheta(:,1).*sinTheta(:,2);
% R(2,2,:) = cosTheta(:,1).*cosTheta(:,3) - sinTheta(:,1).*sinTheta(:,2).*sinTheta(:,3);
% R(2,3,:) = -cosTheta(:,2).*sinTheta(:,1);
% R(3,1,:) = sinTheta(:,1).*sinTheta(:,3) - cosTheta(:,1).*cosTheta(:,3).*sinTheta(:,2);
% R(3,2,:) = cosTheta(:,3).*sinTheta(:,1) + cosTheta(:,1).*sinTheta(:,2).*sinTheta(:,3);
% R(3,3,:) = cosTheta(:,1).*cosTheta(:,2);

%%

% Calculate indices for accessing rotation matrix
i = 3; j = 2; k = 1;

% Find special cases of rotation matrix values that correspond to Euler
% angle singularities.
sy = sqrt(R(i,i,:).*R(i,i,:) + R(j,i,:).*R(j,i,:));
singular = sy < 10 * eps(class(R));

% Calculate Euler angles
a = [-atan2(R(j,i,:), R(i,i,:)),...
     -atan2(-R(k,i,:), sy),...
     -atan2(R(k,j,:), R(k,k,:))];

% Singular matrices need special treatment
numSingular = sum(singular,3);
assert(numSingular <= length(singular));
if numSingular > 0
    a(:,:,singular) = [zeros(1,1,numSingular,'like',R), ...
                       -atan2(-R(k,i,singular), sy(:,:,singular)), ...
                       -atan2(-R(j,k,singular), R(j,j,singular))];
end

