function [L]=ind2logic(varargin)

% function [L]=ind2logic(ind,siz)
% ------------------------------------------------------------------------
% This function converts the linear indicies ind to the logic array L. 
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        ind=varargin{1};
        siz=[max(ind(:)) 1];
    case 2
        ind=varargin{1};
        siz=varargin{2};        
end

%Assume vector if size has one entry
if numel(siz)==1
    siz=[siz 1];
end

%%
% Create logic
L=false(siz); %Initialise as all false
L(ind)=1; %Set to true and the indices ind

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
