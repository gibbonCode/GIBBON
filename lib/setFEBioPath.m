function setFEBioPath(FEBioPathSpec)

%%
filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths

T{1}='#List paths for febio here, multiple can be given e.g. for use on several platforms. The first valid path found is always used.';

%Set paths in cell
switch class(FEBioPathSpec)
    case 'cell'
        for q=1:1:numel(FEBioPathSpec)
            pathNameNow=FEBioPathSpec{q};
            if exist(pathNameNow,'file')==2
                T{q+1}=FEBioPathSpec{q};
            end
        end
    otherwise
        if exist(FEBioPathSpec,'file')==2
            T{2}=FEBioPathSpec;
        end
end

cell2txtfile(fileName,T,0); %Write cell to config file

end


 
%% <-- GIBBON footer text --> 
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
