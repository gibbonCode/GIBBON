function cVoigt=voigtMap(c)

siz2=3*ones(1,2);
siz4=3*ones(1,4);

[linearIndexVoigt,linearIndexFourthOrder]=tensor2voigtMap(c);
secondOrder=0;
if ndims(c)==4 %4d array e.g. 4th order tensor
    if all(size(c)==siz4) %4th order tensor
        switch class(c)
            case 'double'
                cVoigt=zeros(6,6);
            case 'sym'
                cVoigt=sym(zeros(6,6));
        end
    else
        %TO DO add warning here
    end
elseif all(size(c)==siz2) %2nd order tensor
    secondOrder=1;
        switch class(c)
            case 'double'
                cVoigt=zeros(6,1);
            case 'sym'
                cVoigt=sym(zeros(6,1));
        end
else
    %TO DO add warning here
end

cVoigt(linearIndexVoigt)=c(linearIndexFourthOrder);
if secondOrder==1
    cVoigt(4:end)=2.*cVoigt(4:end);
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
