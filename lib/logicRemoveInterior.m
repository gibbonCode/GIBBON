function [Lb]=logicRemoveInterior(L)

% function [Lb]=logicRemoveInterior(L)
% ------------------------------------------------------------------------
% This function removes the interior entries for the input logic L.
% Interior entries are those fully surrounded by neighbours (e.g. in 3D all
% entries which are attached to a top, bottom, left, right, front, and a
% back neighbour).
% Vectors and n-dimensional arrays are supported. 
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% Change log:
% 2018/06/14 Created as alternative to bwmorph3 which requires a special
% toolbox and is from R2018a. 
%------------------------------------------------------------------------

if ~isempty(L)
    if isvector(L) %Handle vector
        h=ones(1,3);
        if ~isrow(L)
            h=h'; %Transpose for column
        end
    else
        nd=ndims(L);
        h=zeros(3*ones(1,nd)); %Initialize h
        ijkMid=2*ones(1,nd); %Middle coordinate
        I=eye(nd,nd); %Identify matrix to offset indices
        ijkFilter=[ijkMid; ijkMid(ones(nd,1),:)+I; ijkMid(ones(nd,1),:)-I]; %Indices for filer
        ind=sub2indn(size(h),ijkFilter);%Linear indices for filer
        h(ind)=1; %Set ones for ND cross shape
    end
    L=L>0;
    LC=convn(double(L),h,'same');
    Lb=L & ~(LC==sum(h(:)));
else
    Lb=[];
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
