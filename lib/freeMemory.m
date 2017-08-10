function [numFreeBytes]=freeMemory

try
    if ispc
        [~,mem_stat] = memory;
        numFreeBytes = mem_stat.PhysicalMemory.Available;
        % numFreeBytes = mem_stat.MaxPossibleArrayBytes; %Alternative
    elseif isunix && ~ismac
        % Output format of running 'free -b | grep Mem'
        %  total       used       free     shared    buffers     cached
        
        %         [~,S] = unix('free -b | grep Mem'); % Excute free command and collect output in strin S
        %         mem_stat = str2double(regexp(S, '[0-9]*', 'match')); % Get the numbers
        %         numFreeBytes = mem_stat(3) + mem_stat(end) ;
        
        [~,numFreeBytesStr]=unix('free -b | awk ''/Mem/{print $3} /Mem/{print $6}''');
        numFreeBytes=sum(str2num(numFreeBytesStr));
    elseif ismac 
        %UNTESTED!
        [~,S] = unix('vm_stat | grep free'); % Excute vm_stat
        mem_stat = strfind(S,' '); % Detect spaces
        numFreeBytes = str2double(S(mem_stat(end):end))*4096; %Take last value and convert pages to bytes
    end
catch 
    error('Could not determine free memory');
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
