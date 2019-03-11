function c=voigtUnMap(cVoigt)


if isvector(cVoigt) %assume that c is a 4th order tensor
    siz_c=[3 3];
    cVoigt(4:end)=(1/2).*cVoigt(4:end); %Undo doubling
else
    siz_c=[3 3 3 3];
end

[linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(zeros(siz_c));
    
switch class(cVoigt)
    case 'double'
        c=zeros(siz_c);
    case 'sym'
        c=sym(zeros(siz_c));
end

c(linearIndexFourthOrder)=cVoigt(linearIndexVoigt);

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
