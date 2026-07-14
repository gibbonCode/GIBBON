function getUndocumented(varargin)

% function getUndocumented(n)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=1;
    case 1
        n=varargin{1};
end

%%
testFolder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'docs');
[testFileList]=getTestFiles('HELP'); 

c=1; 
for q=1:1:numel(testFileList)

    mFileNow=fullfile(testFolder,testFileList{q});

    % fid = fopen(mFileNow);
    % lineNow = fgetl(fid) 

    fileLines = readlines(mFileNow);

    if any(contains(fileLines,'UNDOCUMENTED'))
        if c==n
            open(mFileNow)
            return
        end
        c=c+1;
    end    
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
