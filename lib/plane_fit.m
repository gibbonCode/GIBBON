function [A,B,C]=plane_fit(x,y,z)

% function [A,B,C]=plane_fit(x,y,z)
% ------------------------------------------------------------------------
%   Fit a plane to x,y,z data.
%   [A,B,C]=plane_fit(x,y,z) calculates the coefficients A,B,C that fit the data
%   defined by the vectors x,y,z. 
%
%   Uses command svd
%
%   %EXAMPLE: 
%
%         [x,y]=meshgrid(linspace(0,10,20),linspace(0,10,20));
%         a=1; b=2; c=-2;
%         z=(a*x)+(b*y)+c;
%         x=x(:); y=y(:); z=z(:);
%         z=z+(randn(length(z),1));
%         [A,B,C]=plane_fit(x,y,z); 
%         [X,Y]=meshgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%         Z=(A*X)+(B*Y)+C;
%         plot3(x,y,z,'r.'); hold on; grid on;
%         surf(X,Y,Z,'FaceColor','g'); alpha(0.5);
%         title(['a=',num2str(a), ', A=',num2str(A),', b=',num2str(b),', B=',num2str(B),', c=',num2str(c),', C=',num2str(C)]);
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 14/08/2008
% ------------------------------------------------------------------------

P=[mean(x),mean(y),mean(z)];
[U,S,V]=svd([x-P(1),y-P(2),z-P(3)],0);
N=-1/V(end,end)*V(:,end);
A=N(1); B=N(2); C=-P*N;
 
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
