function [Xc,Yc]=polycentroid(X,Y)

%N.B. 
% Assumes row vectors or matrices whereby each row describes a polygon with
% points appearing in the order defining the polygon

meanX=mean(X,2);
meanY=mean(Y,2);
X=X-meanX*ones(1,size(X,2));
Y=Y-meanY*ones(1,size(Y,2));

A = polyarea(X',Y');
Xc=sum((X(:,1:end-1)+X(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));
Yc=sum((Y(:,1:end-1)+Y(:,2:end)).*(X(:,1:end-1).*Y(:,2:end)-X(:,2:end).*Y(:,1:end-1)),2)./(6*A(:));

Xc=Xc+meanX;
Yc=Yc+meanY;

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
