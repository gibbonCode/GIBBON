function [V]=randInCircle(varargin)

% function [V]=randInCircle(n,R)
% ------------------------------------------------------------------------
% Create random and uniform distribution of n points in a circle with
% radius R. The output array V is an nx2 2D coordinate array. 
%
% 2021/06/28 Created KMM
% ------------------------------------------------------------------------

%%

switch nargin
    case 1
        n=varargin{1};
        R=1;
    case 2
        n=varargin{1};
        R=varargin{2};
    otherwise
        error('Wrong number of input arguments. Define 1 or 2 input arguments'); 
end

%%

r = R.*sqrt(rand(n,1)); %Radii
theta = rand_angle([n,1]); %rand(n,1) * 2*pi; %Angles
V=[r.*cos(theta) r.*sin(theta)]; %Points

