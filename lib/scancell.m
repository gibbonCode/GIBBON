function [IND_found]=scancell(T,targets,found_count)

line_count=1;
IND_found=cell(numel(targets),1);

if isempty(found_count)
    while 1
        for i=1:1:numel(targets)
            if ~isempty(strfind(T{line_count},targets{i}));
                IND_found{i}=[IND_found{i} line_count];
            end
        end
        line_count=line_count+1;
        if line_count>numel(T)
            break
        end
    end
else
    FOUND_count=zeros(size(found_count));
    while 1
        for i=1:1:numel(targets)
            if ~isempty(strfind(T{line_count},targets{i}));
                IND_found{i}=[IND_found{i} line_count];
                FOUND_count(i)=FOUND_count(i)+1;
            end
        end
        line_count=line_count+1;
        if line_count>numel(T)
            break
        end
        if all(FOUND_count==found_count)
            break
        end
    end
end

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
