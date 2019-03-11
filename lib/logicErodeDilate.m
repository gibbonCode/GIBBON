function [Lk]=logicErodeDilate(L,k,erodeDilateOption)

%Normalise structural kernel
k=k./sum(k(:));

%Convolve logic with kernel
M=double(L);
p=size(k)-1;
M_rep=zeros(size(M)+p);
M_rep(p(1)-1+(1:size(M,1)),p(2)-1+(1:size(M,2)),p(3)-1+(1:size(M,3)))=M;
Mk = convn(M_rep,k,'valid');
Mk=Mk./max(Mk(:)); %Scale max to 1

epsMax=max(eps(Mk(:)));

%Erode or dilate logic
switch erodeDilateOption
    case 'erode'
        Lk=(Mk>=(1-epsMax)); %Stayed nearly 1
    case 'dilate'    
        Lk=Mk>=(0+epsMax); %Became higher than 0
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
