function [Vf,indFix]=curvePathOrderFix(V)

%%

%Deriving distance matrix
try
    D=dist(V,V');
catch
    D=distND(V,V);
end

D(eye(size(D))==1)=nan;

%Start loop
indFix=ones(1,size(V,1));
% hw = waitbar(0,'curvePathOrderFix:');
numPoints=size(V,1);
for qIter=2:1:numPoints
%     waitbar(qIter/numPoints);    
    [~,indMin]=gnanmin(D,[],2); 
    indFix(qIter)=indMin(indFix(qIter-1));    
    D(:,indFix)=NaN;                
end
Vf=V(indFix,:);
% close(hw);
 
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
