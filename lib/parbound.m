function Pb=parbound(P,Ps)

% function Pb=parbound(P,Ps)
% ------------------------------------------------------------------------
% Maps the values in the vector P to Pb according to upper and lower limits
% set in the input structure Ps. 
%
% Ps.ub -> Upper bound
% Ps.lb -> Lower bound
% Ps.c -> Centre of the output interval (i.e. initial parameter)
% Ps.f -> Scale factor such that at for the input range
% [-f*(ub-lb),f*(ub-lb)] the width of the output range is Ps.t*(ub-lb)
% Ps.t -> Factor (0-1) of the output range width in relation to the intput
% range width 
%
% 
% %% EXAMPLE
% Ps.f=2; 
% Ps.t=0.9;
% Ps.ub=10;
% Ps.lb=0;
% Ps.c=7;
% 
% w=Ps.ub-Ps.lb;
% P=linspace(Ps.c-Ps.f*w,Ps.c+Ps.f*w,500);
% Pb=parbound(P,Ps);
% hf1=cFigure; hold on; grid on;
% xlabel('P','FontSize'); ylabel('Pb','FontSize'); 
% hf=plot(P,Pb,'b-','LineWidth',5);
% hf=plot(P,(Ps.ub-(0.5*(1-Ps.t)*w))*ones(size(P)),'g-','LineWidth',5);
% hf=plot(P,(Ps.lb+(0.5*(1-Ps.t)*w))*ones(size(P)),'g-','LineWidth',5);
% hf=plot(P,Ps.ub*ones(size(P)),'r-','LineWidth',5);
% hf=plot(P,Ps.lb*ones(size(P)),'r-','LineWidth',5);
% axis equal; axis tight; ylim([Ps.lb Ps.ub]);
% 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 06/02/2012
%------------------------------------------------------------------------

w=abs(Ps.ub-Ps.lb); %Interval width
b=Ps.f*w; %Point at which the output interval width equals t*w
a=tan(Ps.t*0.5*pi);
x=(P-Ps.c)./(b./a);
Pb=(w.*(0.5+(atan(x)./pi)))+Ps.lb;
 
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
