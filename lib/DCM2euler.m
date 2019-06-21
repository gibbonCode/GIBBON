function [a]=DCM2euler(Q)

% function [a]=DCM2euler(Q)
% ------------------------------------------------------------------------
%
%
% Change log:
% 2019/06/21 Simplified code, works for
% symbolic input, assumes non-singular matrices. 
% ------------------------------------------------------------------------

%%

a = [-atan2( Q(2,3,:),  Q(3,3,:)),...
     -atan2(-Q(1,3,:),  sqrt(Q(3,3,:).*Q(3,3,:) + Q(2,3,:).*Q(2,3,:))),...
     -atan2( Q(1,2,:),  Q(1,1,:))];

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
