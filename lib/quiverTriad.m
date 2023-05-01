function [varargout]=quiverTriad(varargin)

% function [h]=quiverTriad(V,R,vecSize,colorOpt,alphaLevel)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        V=varargin{1};
        R=varargin{2};
        vecSize=1;
        colorOpt=[];
        alphaLevel=1;
    case 3
        V=varargin{1};
        R=varargin{2};
        vecSize=varargin{3};
        colorOpt=[];
        alphaLevel=1;
    case 4
        V=varargin{1};
        R=varargin{2};
        vecSize=varargin{3};
        colorOpt=varargin{4};        
        alphaLevel=1;
    case 5        
        V=varargin{1};
        R=varargin{2};
        vecSize=varargin{3};
        colorOpt=varargin{4};
        alphaLevel=varargin{5};
end

%%

%Columns represent vectors

V=V(ones(1,3),:);
[F,V,C]=quiver3Dpatch(V(:,1),V(:,2),V(:,3),R(1,:)',R(2,:)',R(3,:)',eye(3,3),[vecSize vecSize]);

if isempty(colorOpt) %If empty use colormapping
    h=gpatch(F,V,C,'none',alphaLevel);
else %else use specified which could be 'k'
    h=gpatch(F,V,colorOpt,'none',alphaLevel);
end

if nargout==1
    varargout{1}=h;
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
