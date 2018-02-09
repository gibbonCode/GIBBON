function d=divisors(x,arg)

% Determine divisors of input x (which needs to be natural number)
% By default the output provides a vector containing all positive divisors
% but this can be altered depending on arg. 


%Default if no option is given
if exist('arg')==0
    arg='pos';
end

x=abs(x);
X=0:1:x;

d=X(ismember(X,abs(x)./X));
switch arg
    case 'pos'
    case 'neg'
        d=-fliplr(d);
    case 'all'
        d=[-fliplr(d) d];
    otherwise
        error('False input!');
end

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
% Copyright (C) 2018  Kevin Mattheus Moerman
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
