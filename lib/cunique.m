function [varargout]=cunique(varargin)

% function [A_uni,ind1,ind2,Ac]=cunique(A)
%-------------------------------------------------------------------------
%This function is similar to MATLAB's unique function. There are three
%differences: 1) An additional 4th optional output is available providing
%the count, or number of occurances, for each element in the input array.
%2) The 2nd output mathes the size of the first input, 3) The 3rd output
%is reshaped to be the size of the input variable. 
%
% See also: unique
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% Change log: 
% 2018/03/21 Created
% 2019/07/02 Adding varargin handling (pass unique options) including
% 'rows' support. 
% 2020/01/10 Fixed bug in accumarray array size specification
%-------------------------------------------------------------------------

%%

A=varargin{1};

[A_uni,ind1,ind2]=unique(varargin{:});

if any(strcmp(varargin,'rows'))
    typeOpt=2;
else
    typeOpt=1;
end

switch typeOpt
    case 1
        varargout{1}=A_uni;
        varargout{2}=reshape(ind1,size(A_uni));
        varargout{3}=reshape(ind2,size(A));
    case 2
        varargout{1}=A_uni;
        varargout{2}=reshape(ind1,[size(A_uni,1) 1]);
        varargout{3}=reshape(ind2,[size(A,1) 1]);        
end

if nargout==4            
    Ac=accumarray(ind2,1,[length(ind1) 1]);
    Ac=Ac(ind2);
    if typeOpt==1 && size(Ac,1)~=size(A,1)
        Ac=reshape(Ac,size(A));
    end

    varargout{4}=Ac;
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
