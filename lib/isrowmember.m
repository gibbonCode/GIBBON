function logicMember=isrowmember(varargin)

% function logicMember=isrowmember(F1,F2,orderOption)
% ------------------------------------------------------------------------
% 
% ------------------------------------------------------------------------

%% Parse input 

switch nargin    
    case 2
        F1=varargin{1};
        F2=varargin{2};
        orderOption=1;
    case 3
        F1=varargin{1};
        F2=varargin{2};
        orderOption=varargin{3};       
end

%%

%Sort if needed
if orderOption
    %Sort e.g. so that [4 3 2 1]=[1 2 3 4]
    F1=sort(F1,2);
    F2=sort(F2,2);
end

%Create virtual indicesf or each face
maxIndex=max([F1(:);F2(:)]);
sizVirt=maxIndex*ones(1,size(F1,2));
indVirt_F1=sub2indn(sizVirt,F1);
indVirt_F2=sub2indn(sizVirt,F2);

logicMember=ismember(indVirt_F1,indVirt_F2);

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
