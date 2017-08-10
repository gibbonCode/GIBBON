function changeFileExtensions(pathName,extOld,extNew)

%%
if isempty(extOld)
    files = dir(pathName);
    files={files(1:end).name};
else
    files = dir(fullfile(pathName,['*.',extOld]));
    files={files(1:end).name};
end

for q=1:1:numel(files)
    oldNameFull=fullfile(pathName,files{q});
    if ~isdir(oldNameFull)
        [~,oldName,extFile] = fileparts(oldNameFull);
        if isempty(extNew)
            newNameFull=fullfile(pathName,oldName);
        else
            newNameFull=fullfile(pathName,[oldName,'.',extNew]);
        end
        if strcmp(oldNameFull,newNameFull)==0
            if isempty(extOld)
                if strcmp(extFile,extOld);
                    movefileNow(oldNameFull,newNameFull);
                end
            else                
                movefileNow(oldNameFull,newNameFull);
            end
        end
    end
end

end

function movefileNow(oldNameFull,newNameFull)
    movefile(oldNameFull,newNameFull);
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
