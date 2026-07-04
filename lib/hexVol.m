function [varargout]=hexVol(varargin)

% function [VE,logicPositive]=hexVol(E,V,absOpt)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the hexahedral elements specified by the
% element matrix E and the vertices V. The optional input absOpt sets
% wether the volumes are made absolute or if the output may contain
% negative volumes (e.g. for inverted elements). 
% 
% This implementation is based on: 
% https://itk.org/Wiki/images/6/6b/VerdictManual-revA.pdf
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/04/10 KMM: Created 
% 2021/12/10 KMM: Removed reliance on tetrahedron conversion (inaccurate volume
% estimates), and switched to proper hex volume evaluation. 
% 2021/12/10 KMM: Added option to output volume as absolute value or not
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

%Hex vertices
P0=V(E(:,1),:);
P1=V(E(:,2),:);
P2=V(E(:,3),:);
P3=V(E(:,4),:);
P4=V(E(:,5),:);
P5=V(E(:,6),:);
P6=V(E(:,7),:);
P7=V(E(:,8),:);

%Principal axes
X1 = (P1-P0) + (P2-P3) + (P5-P4) + (P6-P7);
X2 = (P3-P0) + (P2-P1) + (P7-P4) + (P6-P5);
X3 = (P4-P0) + (P5-P1) + (P6-P2) + (P7-P3);

% %Compute absolute volumes from 1/64 times the product of the magnitudes of the 3 principal axes
% v=1/64.*sqrt(sum(X1.^2,2)).*sqrt(sum(X2.^2,2)).*sqrt(sum(X3.^2,2));

%Compute volume from determinant of jacobian
a=X1(:,1); b=X1(:,2); c=X1(:,3);
d=X2(:,1); e=X2(:,2); f=X2(:,3);
g=X3(:,1); h=X3(:,2); i=X3(:,3);
v=1/64*((a.*e.*i)+(b.*f.*g)+(c.*d.*h)-(c.*e.*g)-(b.*d.*i)-(a.*f.*h));
   
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
