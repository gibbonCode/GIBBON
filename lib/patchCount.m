function [F_uni,F_count]=patchCount(F)

%Removing double FACES
Fs=sort(F,2); %Sort so faces with same nodes have the same rows
[~,IND_F,IND_F_2]=unique(Fs,'rows');
F_uni=F(IND_F,:);

%Get face counts
numF=size(Fs,1);
numFuni=size(F_uni,1);
logicColourMatrixEntry=sparse(IND_F_2,1:numF,1,numFuni,numF,numF);
F_count=full(sum(logicColourMatrixEntry,2));

 
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
