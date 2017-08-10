function cell2txtfile(fileName,T,skipOpt)

% function cell2txtfile(fileName,T,skipOpt)
% ------------------------------------------------------------------------
%
% This function exports the content in the cell array T to
% the text file fileName. Each entry in the cell array will be a line in
% the txt file. Prior to text file creation the cell is converted to a
% column format. If the input skipOpt=1 cell entries which appear empty
% (after spaces are removed) will be skipped.
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2016/09/09: Updated for GIBBON
% 2016/09/09 added conversion to char and associated warning for
% non-character content. 
%------------------------------------------------------------------------

%%

%Make column
T=T(:);

fid=fopen(fileName,'w');
for q=1:size(T,1)
    l=T{q,:};

    if ~ischar(l) %If this isn't a char then attempt conversion
        l=sprintf('%u',l);
%        warning(['Entry ',num2str(q),' is not a char and was converted to: ',l]);
    end
    
    if skipOpt==0 || ~isempty(deblank(l))
        fprintf(fid,'%s\r\n',l);
    end
    
end
fclose(fid);
 
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
