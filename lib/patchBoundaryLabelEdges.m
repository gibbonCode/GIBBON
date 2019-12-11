function [Eb]=patchBoundaryLabelEdges(F,V,C)

% function [Eb]=patchBoundaryLabelEdges(F,V,C)
% ------------------------------------------------------------------------
% This function outputs the boundary edges for the input patch data F, V,
% with face (e.g. color) labels C. Besides the normal patch boundary the
% boundaries of each labelled set are exported as well. 
% 
% Change log: 
% 2019/07/22 Created
% ------------------------------------------------------------------------

%%

E=patchEdges(F,1);% The patch edges
Es=sort(E,2); % Edges sorted in  column dir so 1 3 2 is the same as 1 2 3
c=unique(C); %Get color labels
logicBoundaryEdges=false(size(E,1),1);
for q=1:1:numel(c)
    Eb_now=sort(patchBoundary(F(C==c(q),:),V),2); %Current sorted boundary edge set    
    logicBoundaryEdges=logicBoundaryEdges | ismember(Es,Eb_now,'rows');
end
Eb=E(logicBoundaryEdges,:);

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
