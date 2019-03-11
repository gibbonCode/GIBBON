function minFig(varargin)

%% Parse input

switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

if isempty(hf)
    hf=gcf; %Take the current figure window handle if the handle is empty
end

%% Turn off warning (JavaFrame will become obsolete, this will lead to error in the future)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

%% Minimize using JaveFrame
try
    pause(0.01);
    jFrame=get(hf,'JavaFrame');
%     set(jFrame,'Minimized',1);
    jFrame.setMinimized(1);
catch exception    
    rethrow(exception)
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
