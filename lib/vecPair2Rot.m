function [R]=vecPair2Rot(varargin)

% function [R]=vecPair2Rot(na,nb)
% ------------------------------------------------------------------------
% This function computes the rotation matrix required to rotate the vector
% na to nb. 
% 
% See also: |vecAngle2Rot|
% 
% 2020/07/02 Created
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1 
        na=varargin{1};        
        nb=[];
    case 2
        na=varargin{1};
        nb=varargin{2};
end

if isempty(nb)
    nb=[0 0 1]; %Use z-vector as default
end

%Force row vectors
if size(na,1)>1
    na=na';
end

if size(nb,1)>1
    nb=nb';
end

%Normalize vectors
na=vecnormalize(na);
nb=vecnormalize(nb);

%%

dn=dot(na,nb); %Dot product of vectors
rotAngle=acos(dn); %Angle between vectors

%Determine rotation axis and formulate rotation matrix
if isapprox(dn,1,eps(1)) % close to na==nb
    R=eye(3,3); %Use identity matrix (no rotation)
else %Significant enough angle
    %Need to compose any orthogonal vector
    if isapprox(dn,-1,eps(1)) %close to n==-nz        
        if ~isapprox(abs(dot(na,[0 1 0])),1,eps(1)) %Try y-axis 
            rotAxis=cross(na,[0 1 0]); %Use x-axis as rotation vector
        else %Find best alternative
            I=eye(3,3); %The identify matrix (contains all coordinate axes)
            dMag=abs(dot(I,na(ones(1,3),:),1)); %Compute dot product with x,y,z
            [~,indMin]=min(dMag); %Get most orthogonal
            rotAxis=cross(na,I(indMin,:)); %Compute orthogonal vector
        end        
    else %Not alligned or anti-alligned n~=nz or n~=-nz
        rotAxis=vecnormalize(cross(na,nb));
    end
    %Compute rotation matrix
    R=vecAngle2Rot(rotAngle,rotAxis);
end
