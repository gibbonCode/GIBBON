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
% hf1=figuremax(fig_color,fig_colordef); hold on; grid on;
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
 
%% <-- GIBBON footer text --> 
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
