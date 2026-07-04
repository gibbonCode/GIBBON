function [V]=bezierCurve(p,n)
% function [V]=bezierCurve(p,n)
% ------------------------------------------------------------------------
% This function uses the control points p to create a Bézier curve. The
% output V consists of n points on the curve. 
% 
% 
% ------------------------------------------------------------------------
%%

% if size(p,2)==2
%     p(:,3)=0;
% end

%%

t=linspace(0,1,n)';
N=size(p,1);
nn=0:1:N-1;
f=factorial(nn); 
s=factorial(N-1)./( f.*flip(f) ); %Sigma
V=( s.* ((1-t).^flip(nn)) .* (t.^nn) )*p; %Output

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
