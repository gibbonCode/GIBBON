function [L_fixed]=triSurfLogicSharpFix(F,L,dirOpt)

% function [L_fixed]=triSurfLogicSharpFix(F,L,dirOpt)
%-------------------------------------------------------------------------
%
% 
%-------------------------------------------------------------------------

%%

if size(F,1)~=size(L,1)
    error('size(F,1)~=size(L,1)');
end

switch dirOpt
    case 1 %Add "inward teeth" to the list
        indInLogic=unique(F(L,:));
        L_fixed=all(ismember(F,indInLogic),2);
    case 2 %Remove sharp stand-alone triangles from list
        indNotInLogic=unique(F(~L,:));
        L_fixed=L & ~all(ismember(F,indNotInLogic),2);                
    case 3 %First 1 then 2
        [L]=triSurfLogicSharpFix(F,L,1);
        [L]=triSurfLogicSharpFix(F,L,2);        
        L_fixed=L;
    case 4 %First 2 then 1
        [L]=triSurfLogicSharpFix(F,L,2);
        [L]=triSurfLogicSharpFix(F,L,1);
        L_fixed=L;
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
