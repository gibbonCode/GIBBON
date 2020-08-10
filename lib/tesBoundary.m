function [indBoundary]=tesBoundary(varargin)

% function [indBoundary]=tesBoundary(F,V)
% ------------------------------------------------------------------------
%
% Change log: 
% 2009
% 2019/04/24 Added varargin support since V can be skipped
% 2019/04/24 Started working on cell (e.g. mixed mesh) support, not
% completed yet 
% 2020/07/09 Added mixed mesh support, works for pentahedra, needs work for
% general case
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        F=varargin{1};
        V=[]; 
    case 2
        F=varargin{1};
        V=varargin{2};
end

if isempty(V)
    if isa(F,'cell')
        numPoints=max(cellfun(@(x) max(x(:)),FE));
    else
        numPoints=max(F(:));
    end
elseif size(V,2)==1 %Assume number of points provided
    numPoints=V;
else %Get number of points from data
    numPoints=size(V,1);
end

%%

if isa(F,'cell')
    indBoundary=cell(size(F));
    for q=1:1:numel(F)
        indBoundary{q}=getBoundary(F{q},numPoints);
    end
else
    [indBoundary]=getBoundary(F,numPoints);
end

end

function [indBoundary]=getBoundary(F,numPoints)
    Fbs=sort(F,2);
    sizVirt=numPoints*ones(1,size(Fbs,2));
    ind_F = sub2indn(sizVirt,Fbs);
    [~,indUni1,~]=unique(Fbs,'rows'); %Get indices for unique faces
    ind_F_uni=ind_F(indUni1,:);
    ind=1:1:size(Fbs,1);
    ind=ind(~ismember(ind,indUni1));
    ind_Fb_cut=ind_F(ind,:);
    L_uni=~ismember(ind_F_uni,ind_Fb_cut);
    indBoundary=indUni1(L_uni,:);
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
