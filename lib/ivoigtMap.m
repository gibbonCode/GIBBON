function c=ivoigtMap(cVoigt)

if ~isvector(cVoigt) %assume that c shall be a 4d array e.g. 4th order tensor
    siz=3*ones(1,4);
elseif isvector(cVoigt) %c shall be a 2D array
    siz=3*ones(1,2); 
    cVoigt(4:end)=cVoigt(4:end)./2;
    secondOrder=1;
else
    %TO DO error catching
end

switch class(cVoigt)
    case 'double'
        c=zeros(siz);
    case 'sym'
        c=sym(zeros(siz));
end

[linearIndexVoigt,linearIndexFull]=tensor2voigtMap(c);
c(linearIndexFull)=cVoigt(linearIndexVoigt);
if secondOrder==1
  [I,J]=ind2sub(siz,linearIndexFull(4:end));
  linearIndexFull2=sub2ind(siz,J,I);
  c(linearIndexFull2)=cVoigt(linearIndexVoigt(4:end));
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
