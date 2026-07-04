function [C]=linspacen(A,B,n)

% function [C]=linspacen(A,B,n)
% ------------------------------------------------------------------------
% This function is a generalization of the linspace function to N
% dimensions. The output C is a matrix of size [size(A) n] such that "it
% goes from A to B in n steps in the last dimention. The input variables A
% and B (scalars, vectors or matrices). For scalar input this function is
% equivalent to linspace.
% The inputs A and B should have the same size.
%
% Change log:
% 2010/07/15 Updated
% 2019/06/29 Fixed bug in relation to numerical precission
% 2019/06/29 Improved error handling
% 2019/06/29 Avoid NaN if n=1, and copy behaviour of linspace for this case
%
%------------------------------------------------------------------------

%%

size_A=size(A); %Store size

if ~all(size(A)==size(B))
    error('A and B should be the same size');
end

if n==1 %Same behaviour as linspace
    C=B;
else
    logicReshape=~isvector(A);
    
    %Make columns
    A=A(:);
    B=B(:);
    
    C=repmat(A,[1,n])+((B-A)./(n-1))*(0:1:n-1);
    
    %Force start and end to match (also avoids numerical precission issues) 
    C(:,1)=A; %Overide start
    C(:,end)=B; %Overide end
    if logicReshape
        C=reshape(C,[size_A n]);
    end
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
