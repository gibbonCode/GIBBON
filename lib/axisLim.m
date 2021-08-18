function axLim=axisLim(V_DEF)

% function axLim=axisLim(V_DEF)
% ------------------------------------------------------------------------
% Computes tight axis limits for the input coordinates V. The coordinates
% may be a 3-dimensional array where the 3rd direction could reflect
% coordinates as a function of time for instance. 
% 
% ------------------------------------------------------------------------
%%

try
    axLim=[min(V_DEF,[],[1 3]); max(V_DEF,[],[1 3])];
catch 
    axLim=[min(min(V_DEF,[],1),[],3); max(max(V_DEF,[],1),[],3)];
end
axLim=axLim(:)';

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
