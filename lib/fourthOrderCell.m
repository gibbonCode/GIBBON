function [CM]=fourthOrderCell(C)

CM=cell(3,3);

for i=1:1:3
    for j=1:1:3
        switch class(C)
            case 'sym'
                C_sub=sym(zeros(3,3));
            otherwise
                C_sub=zeros(3,3);                
        end
        
        for k=1:1:3
            for l=1:1:3                
                C_sub(k,l)=C(i,j,k,l);
            end
        end
        CM{i,j}=C_sub;
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
