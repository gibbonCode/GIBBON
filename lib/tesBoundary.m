function [indBoundary]=tesBoundary(varargin)

% function [indBoundary]=tesBoundary(F)
% ------------------------------------------------------------------------
%
% Change log: 
% 2009
% 2019/04/24 Added varargin support since V can be skipped
% 2019/04/24 Started working on cell (e.g. mixed mesh) support, not
% completed yet 
% 2020/07/09 Added mixed mesh support, works for pentahedra, needs work for
% general case
% 2021/10/08 Simplified
% 2021/10/08 No longer needs vertices as input
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        F=varargin{1};
    case 2
        F=varargin{1};
        warning('Second input (vertices) no longer required. Update code to avoid future error.');
end

%%

if isa(F,'cell')
    indBoundary=cell(size(F));
    for q=1:1:numel(F)
        indBoundary{q}=getBoundary(F{q});
    end
else
    [indBoundary]=getBoundary(F);
end

end

%%

function [indBoundary]=getBoundary(F)

Fs=sort(F,2); %Sort so faces with same nodes have the same rows
[~,~,~,F_count]=cunique(Fs,'rows'); %get indices for unique faces
indBoundary=find(F_count==1);

% Old method
%     Fbs=sort(F,2);
%     sizVirt=numPoints*ones(1,size(Fbs,2));
%     ind_F = sub2indn(sizVirt,Fbs);
%     [~,indUni1,~]=unique(Fbs,'rows'); %Get indices for unique faces
%     ind_F_uni=ind_F(indUni1,:);
%     ind=1:1:size(Fbs,1);
%     ind=ind(~ismember(ind,indUni1));
%     ind_Fb_cut=ind_F(ind,:);
%     L_uni=~ismember(ind_F_uni,ind_Fb_cut);
%     indBoundary=indUni1(L_uni,:);
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
