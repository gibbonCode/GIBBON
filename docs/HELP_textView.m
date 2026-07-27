%% textView
% Below is a demonstration of the features of the |textView| function

%% Syntax
%
%   textView(T);
%   textView(T, optionStruct);
%   textView(T, Name=Value);

%% Description 
% This function displays text content in a scrollable viewer window using a
% monospaced layout with line numbers. Input |T| can be a cell array of text
% lines or a text file path, with optional display settings in a structure.
%
% The following options can be specified in the |optionStruct| or as
% name-value pairs:
%
%   % Name = Default Value
%   BackgroundColor = [44 48 55]/255
%   FontColor = [168 177 190]/255
%   HorizontalAlignment = 'Left'
%   FontName = 'Monospaced'
%   FontWeight = 'Bold'
%   FontSize = 12
%   WordWrap = 'off'
%   Editable = 'off'


%% Examples 
%
% 
%   % Display a cell array of text lines
%   T={'First line';'Second line'; 'A long line to test horizontal scrolling.'};
%   textView(T);
%
%   % View the contents of the textView.m file itself
%   file = which('textView.m');
%   textView(file, 'FontColor', [0,1,0]);
%
%   % Customizing colors and font settings
%   T={'Custom styling, mixing an option structure and name-value pairs.'};
%   optionStruct=struct;
%   optionStruct.BackgroundColor=[0.8 0.6 0.6];
%   optionStruct.FontColor=[0.1 0.3 0.4];
%   optionStruct.WordWrap = 'on';
%   hf=textView(T,optionStruct, FontSize=16, FontName='Courier New');
%


%%
% 
% <<gibbVerySmall.gif>>
% 
% _*GIBBON*_ 
% <www.gibboncode.org>
% 
% _Kevin Mattheus Moerman_, <gibbon.toolbox@gmail.com>
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
