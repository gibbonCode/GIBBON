function [varargout]=runTetView(modelName)

%%

warning('This function will be removed, use the meshView function instead, see HELP_meshView');

%% SETTING TETGEN PATHNAMES

compString=computer; 
switch compString
    case 'PCWIN' %Windows 32-bit
        error('PCWIN 32-bit is not supported. Compile tetGen from the source and alter the code here');
    case 'PCWIN64' %Windows 64-bit
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','win64');
        runNameTetView=fullfile(pathNameTetView,'tetview-win.exe');
    case 'GLNXA64'        
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','lin64');
        runNameTetView=fullfile(pathNameTetView,'tetview-linux');
    case 'MACI64'        
        error('MACI64 is not supported yet. Get TetView online and alter the code here');
    otherwise
        error('Your platform does not seem to be supported. Code your own solution or contact support.')
end

modelName=regexprep(modelName,'\','/');

%% RUN TETVIEW

runString=['"',runNameTetView,'" "',modelName,'" & '];
[runStatus,runCmdHist]=system(runString);

varargout{1}=runStatus;
varargout{2}=runCmdHist;

 
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
