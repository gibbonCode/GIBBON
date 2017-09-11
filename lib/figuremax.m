function [hf]=figuremax(varargin)

% function [fig]=figuremax(Cbg,Cdef)
% ------------------------------------------------------------------------
%
% This function opens a maximized figure window using the color settings
% specified.
% 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 30/09/2014
%------------------------------------------------------------------------

%%

switch nargin
    case 0
        Cbg='w';
        Cdef='white';
    case 1
        Cbg=varargin{1};
        Cdef='white';
    case 2
        Cbg=varargin{1};
        Cdef=varargin{2};
    case 3
        Cbg=varargin{1};
        Cdef=varargin{2};
end

%Open figure with handle
figStruct.Name='GIBBON'; %Figure name
figStruct.Color=Cbg; %Figure background color

%Custom figure properties
figStruct.ColorDef=Cdef; %Setting colordefinitions to black
figStruct.ScreenOffset=0; 

hf=cFigure(figStruct); 

% set(hf,'renderer','OpenGL'); %Default renderer changed, options: painters | zbuffer | OpenGL

%Setting renderer. For RGB and colormap driven There are some bugs for the hardware option, hence changed here
% if ispc
    %opengl software;
    %On UNIX systems, start MATLAB with the command, matlab softwareopengl
% end
 
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
% Copyright (C) 2017  Kevin Mattheus Moerman
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
