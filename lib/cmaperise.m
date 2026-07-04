function [varargout]=cmaperise(varargin)

% function [Cmapped,indMap]=cmaperise(C,cmap,clim)
%-------------------------------------------------------------------------
% This function creates RGB colors for the data C using the colormap cmap
% and the limits clim. 
%
%
% Change log:
% 2019/06/25 Added variable input handling
%-------------------------------------------------------------------------

%% Parse input 

switch nargin
    case 1
        C=varargin{1};
        cmap=[];
        clim=[];
    case 2
        C=varargin{1};
        cmap=varargin{2};
        clim=[];        
    case 3
        C=varargin{1};
        cmap=varargin{2};
        clim=varargin{3};
end

if isempty(cmap)
    cmap=viridis(250);
end

if isempty(clim)
   clim=[min(C(:)) max(C(:))];   
   if diff(clim)<eps(0)
       clim=[min(C(:))-1e-3 max(C(:))+1e-3];
   end
end
    
%% Use colormap to define RGB values

C=C(:); %Color data as column

p=(C-clim(1))./(clim(2)-clim(1));
p(p<0)=0;
p(p>1)=1;

p(isnan(p))=0; %Force nan data to 0

IND_cmap=round((p*(size(cmap,1)-1))+1);
Cmapped=cmap(IND_cmap,:);

%%
varargout{1}=Cmapped;
if nargout==2
    varargout{2}=IND_cmap;
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
