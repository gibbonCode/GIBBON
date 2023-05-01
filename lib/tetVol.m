function [varargout]=tetVol(varargin)

% function [VE,logicPositive]=tetVol(E,V,absOpt)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the tetrahedral elements specified by the
% element matrix E and the vertices V.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/09/25
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

X=V(:,1); Y=V(:,2); Z=V(:,3);
XE=X(E); YE=Y(E); ZE=Z(E);
if size(E,1)==1 %Transpose in this special case
   XE=XE'; YE=YE'; ZE=ZE';
end
A=[XE(:,1) YE(:,1) ZE(:,1)];
B=[XE(:,2) YE(:,2) ZE(:,2)];
C=[XE(:,3) YE(:,3) ZE(:,3)];
D=[XE(:,4) YE(:,4) ZE(:,4)];

v=-dot((A-D),cross((B-D),(C-D),2),2)/6;

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
