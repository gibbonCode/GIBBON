function [B,scrambleIndices]=scramble(varargin)

% function [B,scrambleIndices]=scramble(A,scrambleDim)
% ------------------------------------------------------------------------
% This function scrambles a matrix A allong dimension scrambleDim. The
% scrambling action is based on the randperm function. 
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log:
% 2019/03/28 Created
% 2019/05/03 Fixed bug in relation to scrambling allong a dimension
%------------------------------------------------------------------------

%% Parse input
switch nargin    
    case 1
        A=varargin{1};
        scrambleDim=[];
    case 2
        A=varargin{1};
        scrambleDim=varargin{2};
end

%%

B=A; %Initialize B as A
if isempty(scrambleDim)    
    scrambleIndices=randperm(numel(A));
    B=B(scrambleIndices);
else
    scrambleIndices=randperm(size(A,scrambleDim))';   
    subIND=ind2subn(size(A),1:1:numel(A)); %Subscript index matrix    
    subIND(:,scrambleDim)=scrambleIndices(subIND(:,scrambleDim));
    IND=sub2indn(size(A),subIND); %Linear index matrix
    B(IND)=A(1:end);
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
