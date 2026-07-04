function [P]=triIncenter(F,V)

% function [P]=triIncenter(F,V)
% ------------------------------------------------------------------------
% Computes the triangle incenter for the triangles defined by the face
% array F and the vertex array V. 
% 
% Change log: 
% 2022/04/28 Created
% ------------------------------------------------------------------------

%%

try %Try to use triangulation method
    warning('off','MATLAB:triangulation:PtsNotInTriWarnId');
    TR=triangulation(F,V);
    P=incenter(TR);
    warning('on','MATLAB:triangulation:PtsNotInTriWarnId');
catch %Use custom approach instead
    a=sqrt( sum( (V(F(:,2),:)-V(F(:,3),:)).^2 ,2) );
    b=sqrt( sum( (V(F(:,1),:)-V(F(:,3),:)).^2 ,2) );
    c=sqrt( sum( (V(F(:,1),:)-V(F(:,2),:)).^2 ,2) );

    P=zeros(size(F,1),size(V,2));
    ABC=[a b c];
    sABC=sum(ABC,2);
    for q=1:1:size(V,2)
        X=V(:,q);
        P(:,q)=sum(X(F).*ABC,2)./sABC;
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
