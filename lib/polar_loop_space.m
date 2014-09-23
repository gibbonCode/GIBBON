function [X,Y]=polar_loop_space(x,y,n,THETA_start,interpMethod)

% [X,Y]=polar_loop_space(x,y,n,THETA_start,interp_method)
% ------------------------------------------------------------------------
%This function interpolates n evenly spaced points onto the loop specified
%by x,y using interpMethod  as the interpolation method. The loop should
%be interpolatable in a polar coordinate system, i.e. each angle should
%provide a unique radius. THETA_start indicates the angle for the start of
%the loop. 
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 11/04/2011
%------------------------------------------------------------------------

%% Converting to polar coordinates
x=x(:); y=y(:);
[THETA,R] = cart2pol(x,y);
[THETA,ind_sort]=sort(THETA); %Sorting angles
R=R(ind_sort); %Sorting radii according to angle order
THETA=unwrap(THETA); %unwrapping angles

% Removing double points
[THETA,ia,ic] = unique(THETA); clear ic; %N.B. performance depends on numerical precision
R=R(ia);
% THETA(THETA<0)=THETA(THETA<0)+2*pi;

%Closing loop
THETA(end+1)=THETA(1)+2*pi; R(end+1)=R(1);
% THETA(THETA>pi)=(THETA(THETA>pi)-2*pi);

[x,y] = pol2cart(THETA,R); %Getting corresponding cartesian coordinates


%% Deriving "curve distance" vector
D=pathLength([x y]);

%% Deriving distance at THETA_start and THETA_end
% THETA_start=rem(THETA_start,pi);
% THETA_end=rem(THETA_end,pi);

% min(THETA)./pi
% max(THETA)./pi
D_start=interp1(THETA,D,THETA_start,interpMethod);
% D_end=interp1(THETA,D,THETA_end,interp_method);
% if D_start==D_end
%    D_end=D_end+max(D(:)); 
% end
% D_start
% D_end
%% Interpolating "geodesic" points
D_geo=linspace(D_start,D_start+max(D(:)),n+1)'; %n+1 since start=end
D_geo(D_geo>max(D))=D_geo(D_geo>max(D))-max(D);

THETA_geo=interp1(D,THETA,D_geo,interpMethod); %Interpolating angles
R_geo=interp1(D,R,D_geo,interpMethod); %Interpolating radii
[X,Y] = pol2cart(THETA_geo,R_geo); %Getting Cartesian coordinates
X=X(1:end-1); Y=Y(1:end-1); %Removing double end points (=start point)

end

