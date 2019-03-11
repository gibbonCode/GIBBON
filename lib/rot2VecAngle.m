function [theta,w]=rot2VecAngle(R)

% function [theta,w]=rot2VecAngle(R)
% ------------------------------------------------------------------------
%
% 2015
% 2018/11/08 Fixed bug in relation to R close to eye(3,3) causing a zero
% angle and an all NaN rotation vector. 
% ------------------------------------------------------------------------

%%
if any(isnan(R(:)))
    error('Invalid rotation matrix. Matrix contains NaNs');
end

theta=acos(0.5*(trace(R)-1));
theta=real(theta);
w=(1./(2*sin(theta)))*([R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)]);
w=vecnormalize(w);

%Treat devide by zero
if any(isnan(w))
   w=[1 0 0]';
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
