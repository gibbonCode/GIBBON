function [varargout]=batman(n)

% function [varargout]=batman(n)
% ------------------------------------------------------------------------
% The |batman| function implements a particular version of the so called
% batmat-equation, a curve defining the batman logo. The input for this
% function is the number of desired points n. The user may request a sigle
% nx2 output array or two nx1 arrays (x and y coordinates). 
%
% This is a MATLAB implementation of the parameterised Batman equation
% presented by Jerome White (http://www.talljerome.com/ @talljerome,
% https://youtu.be/oaIsCJw0QG8), in particular the form presented here: 
% https://www.desmos.com/calculator/ajnzwedvql
% Modification: The batman is scaled to be 2 in width. 
% 
% 2020/05/05
% ------------------------------------------------------------------------
%%


t=linspace(0,8,1000)';

t1=abs(t-1);
t2=abs(t-2);
t3=abs(t-3);
t4=abs(t-4);
t5=abs(t-5);
t7=abs(t-7);
t8=abs(t-8);

x=(0.3*t)+(0.2*t1)+(2.2*t2)-(2.7*t3)-(3*t5)+(3*t7)...
    +5*sin((pi/4)*(t3-t4+1))+((5/4)*((t4-t5-1)).^3)...
    -(5.3*cos((pi/2+asin(47/53)).*(1/2*(t7-t8-1))))+2.8;
y= 3/2*(t1-t2)+29/4*(t5-t4)+7/16*(t2-t3-1).^4 ...
    +4.5*sin(pi/4*(t3-t4-1))-(3*sqrt(2)/5)*abs(t5-t7).^(5/2)+ ...
    6.4*sin((pi/2+asin(47/53)).*(1/2*(t7-t8+1))+(asin(56/64)))+4.95;

V=[x y;flipud(-x(2:end-1)) flipud(y(2:end-1))]/22;
Vn=evenlySampleCurve(V,n,'pchip',1);

switch nargout
    case 1
        varargout{1}=Vn;
    case 2
        varargout{1}=Vn(:,1);
        varargout{2}=Vn(:,2);
end

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
