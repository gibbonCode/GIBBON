function [varargout]=triPolyDualRefine(varargin)

% function [Ft,Vt,Ct,indIni]=triPolyDualRefine(F,V,fixBoundaryOpt)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/05/01 Updated for GIBBON
% 2015/05/01 Added varargin/varargout
% 2018/04/06 Added fixBoundaryOpt to be compatible with update in
% patch_dual function
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 2
        F=varargin{1}; %Faces
        V=varargin{2}; %Vertices
        fixBoundaryOpt=1;
    case 3        
        F=varargin{1}; %Faces
        V=varargin{2}; %Vertices
        fixBoundaryOpt=varargin{3}; %Vertices
end

%%

[Vd,~,Fds]=patch_dual(V,F,fixBoundaryOpt);

Cd=(1:1:size(Fds,1))';
Cds=Cd(:,ones(size(Fds,2),1));

Fds_t=Fds'; 
Cds_t=Cds';

Vt=[Vd;V];

Q1=Fds_t(:);
L1=Q1>0; 
Q1=Q1(L1);
Ct=Cds_t(:);
Ct=Ct(L1);

Q2=Fds(:,[2:size(Fds,2) 1])';
Q2=Q2(:);
L2=Q2>0; 
Q2=Q2(L2);

q=1:size(Fds,1);
Q3=q(ones(size(Fds,2),1),:)+size(Vd,1);
Q3=Q3(:);
A=Q3;
Q3=Q3(L1);
Ft=[Q1 Q2 Q3];

indIni=(size(Vd,1)+1):(size(Vd,1)+size(V,1));

%% Collect output

varargout{1}=Ft;
varargout{2}=Vt;
varargout{3}=Ct;
varargout{4}=indIni;
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
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
