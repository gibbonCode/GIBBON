function [F,V]=removeNotIndexed(F,V)

% [F,V]=removeNotIndexed(F,V)
% -----------------------------------------------------------------------
% This function takes in a set of indices F (a vector or matrix containing
% indices) into the rows of V and crops V such that only the rows present
% in F are kept. Then the indices in F are fixed such that it contains
% indices for the newly cropped V.
% This function is not specific to patch data but helps in interpretation.
% The names F and V are derived from a patch data type faces (F) and
% vertices array (V). Where the rows in F denotes faces defined by row
% entries in V. Hence this function can be used to remove unused vertices
% from V and to subsequently correct the indices in F. Similarly when
% interpretet in the realm of tesselations F may also define an element or
% edge array.
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/14
%------------------------------------------------------------------------

%%

indUsed=unique(F(:)); %The list of used row indices

%If the number of used entries is not equal to the number of rows in V then
%unused points exist.
if numel(indUsed)<size(V,1)
    %Prepare for index fixing
    indFix=nan(size(V,1),1);
    indFix(indUsed)=1:numel(indUsed);
    
    V=V(indUsed,:); %Crop V to keep only used rows
    F=indFix(F); %Fix indives
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
