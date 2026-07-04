function [P]=tesSmooth(TES,V,IND_V,cPar)

% function [P]=tesSmooth(TES,V,IND_V,cPar)
% ------------------------------------------------------------------------
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/02
%------------------------------------------------------------------------

%% 

%Get/set method
if isfield(cPar,'Method')
    smoothMethod=cPar.Method;
else
    smoothMethod='LAP'; %DEFAULT
end

%Smooth
switch smoothMethod
    case 'LAP'
        [P]=tesSmooth_LAP(TES,V,IND_V,cPar);
    case 'HC'
        [P]=tesSmooth_HC(TES,V,IND_V,cPar);
    otherwise
        error('Invalid smooth method specified');
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
