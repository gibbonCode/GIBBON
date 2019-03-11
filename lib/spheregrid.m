function [X,Y,Z]=spheregrid(xc,yc,zc,n,rc)

% function [X,Y,Z]=spheregrid(xc,yc,zc,n,rc)
% ------------------------------------------------------------------------
% This function generates X,Y,Z coordinates of points inside a circle for
% the 2D case or inside a sphere for the 3D case with a radius of rc. The
% points are distrubuted according to n which is the spacing between points
% in radians.
%
% x,y,z are the x,y,z coordinates of the desired centre point
% n is the resolution of the rotation in radians
% rc is the radius of the required circle/sphere
%
% [X,Y]=spheregrid(x,y,z,revolve_res,r) % will generate points in a circle
% [X,Y,Z]=spheregrid(xc,yc,zc,n,rc) % will generate points in a sphere
%
% %EXAMPLE
% clear all; close all; clc;
% x=2; y=2; z=0; r=1; revolve_res=0.1;
% %2D: points in a circle
% [X,Y]=spheregrid(x,y,z,revolve_res,r);
% figure; fig=gcf; clf(fig); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% subplot(1,2,1);
% plot(x,y,'ro'); hold on; axis equal;
% plot(X,Y,'k.'); hold on; axis equal;
% xlabel('x'); ylabel('y');
% title('2D: points in a circle');
% %3D: points in a sphere
% [X,Y,Z]=spheregrid(x,y,z,revolve_res,r);
% plot3(x,y,z,'ro'); hold on; axis equal;
% [X_sp,Y_sp,Z_sp] = sphere(100); X_sp=(r*X_sp)+x; Y_sp=(r*Y_sp)+y; Z_sp=(r*Z_sp)+z;
% subplot(1,2,2);
% plot3(X,Y,Z,'k.'); hold on; axis equal; grid on;
% surf(X_sp,Y_sp,Z_sp,'EdgeColor','none','FaceColor','g','FaceAlpha',0.5);
% axis equal;hold on;lightangle(80, -40);lightangle(-90, 60);
% title('3D: points in a sphere');

% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 05/05/2009
% ------------------------------------------------------------------------

if nargout==3
    rot_angle=pi;
else
    rot_angle=2*pi;
end

if nargin == 1
    y = x;
end

%Revolve around Z axis
p=floor(rc/n)+1;
x=linspace(0,rc,p);
y=0;
z=0;

[theta, phi, r]=cart2sph(x,zeros(size(x)),z);
THETA=[];
PHI=[];
R=[];
for i= 1: length(r)
    dist_axis=x(i);
    theta_inc = 2*asin((n/2)/dist_axis);
    if isnan(theta_inc)==1
        no_steps=1;
        theta_step=0;
    else
        no_steps=floor(rot_angle/theta_inc)+1;
        theta_step=linspace(0,rot_angle,no_steps)';
        theta_step=theta_step(1:end-1);
    end
    THETA=[THETA; theta_step];
    PHI=[PHI; (phi(i)*ones(1,length(theta_step)))'];
    R=[R; (r(i)*ones(1,length(theta_step)))'];
end
[X Y Z]=sph2cart(THETA,PHI,R);

if nargout==3
    %Revolve around X axis
    [theta, phi, r]=cart2sph(Y,Z,X);
    THETA=[];
    PHI=[];
    R=[];
    rot_angle=2*pi;
    for i= 1: length(r)
        dist_axis=Y(i);
        theta_inc = 2*asin((n/2)/dist_axis);
        if isnan(theta_inc)==1
            no_steps=1;
            theta_step=0;
        else
            no_steps=floor(rot_angle/theta_inc)+1;
            theta_step=linspace(0,rot_angle,no_steps)';
            theta_step=theta_step(1:end-1);
        end
        THETA=[THETA; theta_step];
        PHI=[PHI; (phi(i)*ones(1,length(theta_step)))'];
        R=[R; (r(i)*ones(1,length(theta_step)))'];
    end
    [Y Z X]=sph2cart(THETA,PHI,R);
end

X=X+xc;
Y=Y+yc;
Z=Z+zc;
 
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
