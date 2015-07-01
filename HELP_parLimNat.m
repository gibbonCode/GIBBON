%% parLimNat
% Below is a demonstration of the features of the |parLimNat| function

%% Syntax
% |[xx,S]=parLimNat(xx_c,[xx_min xx_max],x);|

%% Description
% The |parLimNat| function can be used to constrain parameters from
% [-inf,inf] to the range [xx_min xx_max] with the values xx_c at its
% centre. The constraining can be used in combination with optimization
% routinges that do not naturally handle parameter constraints. 

%% Examples

clear; close all; clc;

%%
% PLOT SETTINGS
fontSize=25;
fontSize2=15;
markerSize=45;
lineWidth1=2;
lineWidth2=1;

%% Example: Constraining parameters (normal centre)
xx_c=5;
xx_min=4;
xx_max=10;
x=linspace(xx_c-10,xx_c+10,250);
[xx,S]=parLimNat(xx_c,[xx_min xx_max],x);

hf1=cFigure; hold on; grid on;
title('Constraining using parLimNat','FontSize',fontSize);
xlabel('"free" x','FontSize',fontSize); ylabel('constrained x','FontSize',fontSize); 

hf=plot(x,xx,'r-','LineWidth',lineWidth1);
hf=plot(xx_min,xx_min,'k.','markerSize',markerSize);
hf=plot(xx_max,xx_max,'k.','markerSize',markerSize);
plot(xx_c,xx_c,'k.','markerSize',markerSize);

text(xx_c+0.5,xx_c,'xx_c = centre','Interpreter','none','FontSize',fontSize2);
text(xx_max,xx_max-0.5,'xx_max = upper bound','Interpreter','none','FontSize',fontSize2);
text(xx_min-2.5,xx_min+0.5,'xx_min = lower bound','Interpreter','none','FontSize',fontSize2);

axis tight; axis equal;
set(gca,'FontSize',fontSize);
drawnow;

%% Example: Constraining parameters (out of centre centre)
xx_c=2;
xx_min=0;
xx_max=10;
x=linspace(xx_c-10,xx_c+15,100);
[xx,S]=parLimNat(xx_c,[xx_min xx_max],x);

hf1=cFigure; hold on; grid on;
title('Constraining using parLimNat','FontSize',fontSize);
xlabel('"free" x','FontSize',fontSize); ylabel('constrained x','FontSize',fontSize); 

hf=plot(x,xx,'r-','LineWidth',lineWidth1);
hf=plot(xx_min,xx_min,'k.','markerSize',markerSize);
hf=plot(xx_max,xx_max,'k.','markerSize',markerSize);
plot(xx_c,xx_c,'k.','markerSize',markerSize);

text(xx_c+0.5,xx_c,'xx_c = centre','Interpreter','none','FontSize',fontSize2);
text(xx_max,xx_max-0.5,'xx_max = upper bound','Interpreter','none','FontSize',fontSize2);
text(xx_min-2.5,xx_min+0.5,'xx_min = lower bound','Interpreter','none','FontSize',fontSize2);

axis tight; axis equal;
set(gca,'FontSize',fontSize);
drawnow;

%% 
%
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>