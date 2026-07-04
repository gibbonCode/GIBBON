function [y]=gnansum(varargin)
% function [y]=gnansum(A,[],dimDir)
%------------------------------------------------------------------------
% this function is a replacement for nansum which is part of the
% statistics and machine learning toolbox. 
% The use of sum with the omitnan flag is used or if this does not work
% (for old MATLAB versions) a custom implementation is used. 
%
% Change log:
% 2020/01/09 Created
%------------------------------------------------------------------------

%%
if nargin==1
    varargin{2}=[];
end

try 
    %Use 'omitnan' flag    
    if nargin==1
        varargin{2}=[];
    end
    y = sum(varargin{:},'omitnan');    
catch
    A=varargin{1}; 
    A(isnan(A))=0;
    varargin{1}=A;
    y = sum(varargin{:});    
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
