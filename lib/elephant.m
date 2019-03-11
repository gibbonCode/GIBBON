function [F,V]=elephant

% function [F,V]=elephant
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for an
% elephant model.
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------
%%

% filePath=mfilename('fullpath');
% toolboxPath=fileparts(fileparts(filePath));
% offPath=fullfile(toolboxPath,'data','OFF');
% fileName=fullfile(offPath,'elephant-50kv.off');
% [F,V] = import_off(fileName);
% 
% [F,V,~,~]=triSurfRemoveThreeConnect(F,V,[]);
% 
% [R,~]=euler2DCM([1/3*pi 0 0]);
% V=(R*V')';
% [R,~]=euler2DCM([0 1/36*pi 0]);
% V=(R*V')';

[F,V]=graphicsModels(7);

 
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
