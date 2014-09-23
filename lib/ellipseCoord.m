function [V]=ellipseCoord(A,t)

% function [V]=ellipseCoord(A,t)
% ------------------------------------------------------------------------
% Calculates ellipse coordiantes for the angles in t based on the vector A
% which defines the centre coordinates, the radii and the angle
% respectively. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2013/24/09
%------------------------------------------------------------------------

x0=A(1);
y0=A(2);
x=A(3).*cos(t);
y=A(4).*sin(t);
V=[x(:) y(:) zeros(size(x(:)))];
[R,~]=euler2DCM([0 0 -A(5)]);
V=V*R;
V(:,1)=V(:,1)+x0;
V(:,2)=V(:,2)+y0;