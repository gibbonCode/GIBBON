function [Fn,Vn]=subTriCentre(F,V,L)

% function function [Fn,Vn]=subTriCentre(F,V,L)
% ------------------------------------------------------------------------
% This function splits the faces defined by L up into three by introducing
% a central node. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/03
%------------------------------------------------------------------------

numPoints=size(V,1);
Fs=F(L,:); %Faces to split

%Get face centres
X=V(:,1); Y=V(:,2); Z=V(:,3);
if nnz(L)==1 %treat different behaviour for since face problem
    XF=mean(X(Fs),1);
    YF=mean(Y(Fs),1);
    ZF=mean(Z(Fs),1);
else
    XF=mean(X(Fs),2);
    YF=mean(Y(Fs),2);
    ZF=mean(Z(Fs),2);
end
V_add=[XF(:) YF(:) ZF(:)];
numPointsAdd=size(V_add,1);

Vn=[V;V_add];
indNew=(numPoints+1:numPoints+numPointsAdd)';
F_add=[Fs(:,2) indNew Fs(:,1);...
       Fs(:,3) indNew Fs(:,2);...
       Fs(:,1) indNew Fs(:,3)];

Fn=[F(~L,:);F_add];
 
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
