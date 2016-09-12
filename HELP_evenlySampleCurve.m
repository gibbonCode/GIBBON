%% evenlySampleCurve
% Below is a basic demonstration of the features of the |evenlySampleCurve| function.
%% ind2patch
% Below is a demonstration of the features of the |ind2patch| function

%%
close all; clc; clear;

%% Syntax
% |[Vg] = evenlySampleCurve(V,n,interpPar,closeLoopOpt);|

%% Description
% The |evenlySampleCurve| function samples a curve evenly in n points. The
% curve is parameterized using curve distance and can be closed loop if
% closeLoopOpt==1 (default=0). The resampling is performed using
% interpolation based on the method specified by interpPar. 
% Available methods are those associated with interp1 i.e.: 'linear',
% 'nearest', 'next', 'previous', 'spline', 'pchip' (default), 'cubic'.
% Alternatively interpPar my be set as a scalar in the range 0-1 to use the
% csaps method for cubic spline based smoothening.
%
% See also: |interp1|, |csaps|

%% Examples

%%
% Plot settings
markerSize=15;
lineWidth=2;

%% EXAMPLE USING NORMAL INTERPOLATION

%Simulating the case of an unevenly sampled loop curve
ns=150;
t=sort(linspace(0,2*pi,ns)+pi/10*rand(1,ns));
t=unique(t); %removing double points
t=t(t<2*pi);%Removing 2*pi points since they are the same as the 0 point
r=3+2.*sin(5*t);
[x,y] = pol2cart(t,r);
z=y;
V=[x(:) y(:) z(:)];

interpMethod='pchip';
closeLoopOpt=1;
n=200;
[Vg]=evenlySampleCurve(V,n,interpMethod,closeLoopOpt);

hf1=cFigure;
subplot(1,2,1); hold on;
title('Unevenly sampled');
plot3(V(:,1),V(:,2),V(:,3),'r.-','MarkerSize',markerSize);
drawnow; view(3); grid on; axis equal; axis tight; 
subplot(1,2,2); hold on;
title('Evenly sampled allong curve');
plot3(Vg(:,1),Vg(:,2),Vg(:,3),'g.-','MarkerSize',markerSize);
plot3(Vg(1,1),Vg(1,2),Vg(1,3),'r.','MarkerSize',2*markerSize,'lineWidth',lineWidth);
plot3(Vg(end,1),Vg(end,2),Vg(end,3),'b.','MarkerSize',2*markerSize,'lineWidth',lineWidth);
drawnow; view(3); grid on; axis equal; axis tight; 

%% EXAMPLE USING CSAPS

%Adding noise
V=V+0.2.*randn(size(V));

interpMethod=0.7;
closeLoopOpt=1;
[Vg]=evenlySampleCurve(V,n,interpMethod,closeLoopOpt);

hf1=cFigure;
subplot(1,2,1); hold on;
title('Unevenly sampled');
plot3(V(:,1),V(:,2),V(:,3),'r.-','MarkerSize',markerSize);
drawnow; view(3); grid on; axis equal; axis tight; 
subplot(1,2,2); hold on;
title('Evenly sampled allong curve and smoothened');
plot3(Vg(:,1),Vg(:,2),Vg(:,3),'g.-','MarkerSize',markerSize,'lineWidth',lineWidth);
plot3(Vg(1,1),Vg(1,2),Vg(1,3),'r.','MarkerSize',2*markerSize,'lineWidth',lineWidth);
plot3(Vg(end,1),Vg(end,2),Vg(end,3),'b.','MarkerSize',2*markerSize,'lineWidth',lineWidth);
drawnow; view(3); grid on; axis equal; axis tight; 

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
