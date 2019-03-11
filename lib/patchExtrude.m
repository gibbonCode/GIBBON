function [Et,Vt]=patchExtrude(F,V,zRange)

numSteps=numel(zRange);

%Deriving coordinates
Vt=repmat(V,numSteps,1);
Z_add=ones(size(V,1),1)*zRange; 
Z_add=Z_add(:);
Vt(:,3)=Vt(:,3)+Z_add;

%Replicated faces matrix
F_rep=repmat(F,numSteps-1,1);

%Fix indices since points are copied also
indFix=0:(numSteps-2); 
indFix=indFix(ones(1,size(F,1)),:);
indFix=indFix(:);
indFix=indFix(:,ones(1,size(F_rep,2)));

%Create element matrix
Et=[F_rep+(size(V,1)*indFix) F_rep+(size(V,1)*(indFix+1))];
 
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
