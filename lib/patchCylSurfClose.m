function [varargout]=patchCylSurfClose(varargin)

% function [F,V,C]=patchCylSurfClose(X,Y,Z,C)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 3
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        C=[];
    case 4
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        C=varargin{4};
end

if isempty(C)
    C=ones(size(Z));
end

%%

[F,V,C]=surf2patch(X,Y,Z,C);
I=[(1:size(Z,1)-1)' (1:size(Z,1)-1)' (2:size(Z,1))' (2:size(Z,1))'];
J=[size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1)];
F_sub=sub2ind(size(Z),I,J);
F=[F;F_sub];

%% Collect output

varargout{1}=F;
varargout{2}=V;
if nargout==3
    C=vertexToFaceMeasure(F,C);
    C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5;
    varargout{3}=C;
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
