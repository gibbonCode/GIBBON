function gdoc(varargin)

% function gdoc(functionName)
%------------------------------------------------------------------------
% View XML files using a brower embedded in a figure window or within
% MATLAB brower.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/05/13 %Updated, renamed from febView to xmlView
% 2017/08/18 %Created varargin with viewer option. Create temp file to
% allow viewing of domNode.
%------------------------------------------------------------------------
%%

switch nargin
    case 1
        topic=varargin{1};        
    case 0
        topic='GIBBON';
end

documentationPath=fullfile(fileparts(fileparts(mfilename('fullpath'))),'docs','html');

helpCheckName=fullfile(documentationPath,['HELP_',topic,'.html']);
demoCheckName=fullfile(documentationPath,['DEMO_',topic,'.html']);

if exist(helpCheckName,'file')
    web(helpCheckName, '-helpbrowser');
elseif exist(demoCheckName,'file')
    web(demoCheckName, '-helpbrowser');
else
    disp('No matching HELP or DEMO files found. Searching all documentation instead');
    docsearch(topic);
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
