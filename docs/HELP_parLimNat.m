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
xx_min=0;
xx_max=10;
x=linspace(xx_c-10,xx_c+10,1000);
[xx,S]=parLimNat(xx_c,[xx_min xx_max],x);

cFigure; hold on; grid on;
title('Constraining using parLimNat','FontSize',fontSize);
xlabel('"free" x','FontSize',fontSize); ylabel('constrained x','FontSize',fontSize); 

plot(x,xx,'r-','LineWidth',lineWidth1);
plot(xx_min,xx_min,'k.','markerSize',markerSize);
plot(xx_max,xx_max,'k.','markerSize',markerSize);
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

cFigure; hold on; grid on;
title('Constraining using parLimNat','FontSize',fontSize);
xlabel('"free" x','FontSize',fontSize); ylabel('constrained x','FontSize',fontSize); 

plot(x,xx,'r-','LineWidth',lineWidth1);
plot(xx_min,xx_min,'k.','markerSize',markerSize);
plot(xx_max,xx_max,'k.','markerSize',markerSize);
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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
