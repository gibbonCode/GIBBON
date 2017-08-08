function [X,Y,Z]=polyrevolve(x,z,n)

% function [X,Y,Z]=polyrevolve(x,z,n)
% ------------------------------------------------------------------------
% This function revolves a 2D polygon around the Z-axis.
%
% x is a vector with x coordinates
% z is a vector with z coordinates
% n is the resolution of the rotation in radians
%
% The vectors x and z are converted to spherical coordinates (rho,r). The
% points are then 'copied' around the z-axis. The space between each point
% is defined by n.
%
% %EXAMPLE
% clear all; close all; clc;
% x=0:0.2:2*pi;
% z=sin(x+(pi/2));
% revolve_res=0.2;
% [X,Y,Z]=polyrevolve(x,z,revolve_res);
%
% grid_res=0.05;
% [XI,YI] = meshgrid(min(min(X)):grid_res:max(max(X)),min(min(Y)):grid_res:max(max(Y)));
% ZI = griddata(X,Y,Z,XI,YI,'cubic');
% figure; fig=gcf; clf(fig); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% plot(x,z,'r-'); hold on; axis equal;
% xlabel('x (mm)'); ylabel('z (mm)');
% title('The polygon');
% xlabel('x (mm)'); ylabel('z (mm)');
% figure; fig=gcf; clf(fig); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% plot3(X,Y,Z,'k.'); hold on; axis equal;
% title('Data "copied" around z-axis using POLYREVOLVE');
% xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
% figure; fig=gcf; clf(fig); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% surf(XI,YI,ZI,'EdgeColor','none'); hold on; axis equal;
% shading interp; material shiny; lightangle(45,30);
% title('Surface fit of the revolved data points')
% xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
%
%%EXAMPLE 2
% clear all; close all; clc;
% x=[40 65];
% z=0;
% no_markers=6;
% revolve_res=2*x(1)*sin(pi/no_markers);
% [X,Y,Z]=polyrevolve(x,z,revolve_res);
% Xr=round(X);
% Yr=round(Y);
% Zr=round(Z);
% figure; fig=gcf; clf(fig); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% h=plot(X,Y,'r.'); hold on; axis equal; axis([min(X)-1 max(X)+1 min(Y)-1 max(Y)+1]);
% set(h,'MarkerSize',20);
% h=plot(Xr,Yr,'b.'); hold on; axis equal; axis([min(X)-1 max(X)+1 min(Y)-1 max(Y)+1]);
% set(h,'MarkerSize',20);
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 14/08/2008
% ------------------------------------------------------------------------

[theta, rho, r]=cart2sph(x,zeros(size(x)),z);
THETA=[];
RHO=[];
R=[];
for i= 1: length(r)
    if n>(2*x(i))
        %The distance between points on a circle can not exceed the
        %diameter. These points are not revolved.
        disp(['Warning: n>(2*r(i)), n exceeds maximum value for point no ',num2str(i) ,', this point will not be revolved! '])
        if x(i)==0 %If distance to rotational axis keep point
            no_steps=1;
            theta_step=0;
            THETA=[THETA; theta_step'];
            RHO=[RHO; (rho(i)*ones(1,length(theta_step)))'];
            R=[R; (r(i)*ones(1,length(theta_step)))'];
        end
    else
        theta_inc = real(2*asin((n/2)/x(i))); %Angular spacing between points in radians
        no_steps=round((2*pi)/theta_inc)+1; 
        theta_step=linspace(0,(2*pi),no_steps); 
        theta_step=theta_step(1:end-1);
        THETA=[THETA; theta_step'];
        RHO=[RHO; (rho(i)*ones(1,length(theta_step)))'];
        R=[R; (r(i)*ones(1,length(theta_step)))'];
    end
end

[X Y Z]=sph2cart(THETA,RHO,R);
 
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
