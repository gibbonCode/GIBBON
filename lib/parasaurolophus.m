function [F,V]=parasaurolophus

% [F,V]=parasaurolophus
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% triangulated parasaurolophus dinosaur model. 
%
% This adjusted model consists of 892 triangular faces and 448 vertices. 
% 
% The model was constructed based on the model given here: 
% https://www.rocq.inria.fr/gamma/gamma/download/affichage.php?dir=DINOSAUR/&name=Parasaurolophus
% 
%
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/04/25 Updated for GIBBON
%------------------------------------------------------------------------
%%

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');
fileName=fullfile(pathName,'parasaurolophus.mat');
D=load(fileName);
F=D.F;
V=D.V;
 
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
