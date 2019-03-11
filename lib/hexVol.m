function [varargout]=hexVol(E,V)

% function VE=hexVol(E,V)
% ------------------------------------------------------------------------
% Calculates the volume (VE) of the hexahedral elements specified by the
% element matrix E and the vertices V.
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2018/04/10
%-------------------------------------------------------------------------

%%

[Et,Vt]=hex2tet(E,V,[],5);
[VE,logicPositive]=tetVol(Et,Vt);

if size(E,1)==1
    VE=sum(VE);
    logicPositive=all(logicPositive);
else   
    VE=sum(reshape(VE,[5 size(E,1)])',2);    
    logicPositive=all(reshape(logicPositive,[5 size(E,1)])',2);
end

%%

varargout{1}=abs(VE);
varargout{2}=logicPositive;
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
