function [varargout]=pentaVol(varargin)

% function [VE,logicPositive]=pentaVol(E,V,absOpt)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the pentahedral elements specified by the
% element matrix E and the vertices V. The optional input absOpt sets
% wether the volumes are made absolute or if the output may contain
% negative volumes (e.g. for inverted elements). 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2022/08/03 KMM: Created 
% To do: Remove reliance on tetrahedron conversion (inaccurate volume)
%-------------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};
        absOpt=1;
    case 3        
        E=varargin{1};
        V=varargin{2};
        absOpt=varargin{3};
end

%% Compute element volumes 

v_tet=tetVol([E(:,[1 2 3 4]); E(:,[2 3 4 5]); E(:,[3 4 5 6]);],V); %Volume of tetrahedra
v=sum(reshape(v_tet,[3,size(E,1)]),1)'; %Volume of pentahedra

%% Collect output
if absOpt==1
    varargout{1}=abs(v);
else
    varargout{1}=v;
end

if nargout>1
    varargout{2}=v>0;
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
