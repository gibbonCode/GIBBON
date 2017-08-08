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
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
