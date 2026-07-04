function cell2txtfile(varargin)

% function cell2txtfile(fileName,T,skipOpt,checkChar)
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
% Change log: 
% 2016/09/09: Updated for GIBBON
% 2016/09/09: Added conversion to char and associated warning for non-character content. 
% 2018/02/14: Added whole cell conversion when skipOp==0
% 2018/02/14: Added varargin and defaults
%
% To do: 
% * Proper handling of conversion of non-character entries (e.g. have user
% specify conversion format to use for numerical data). 
%------------------------------------------------------------------------

%% parse input

switch nargin
    case 2
        fileName=varargin{1};
        T=varargin{2};
        skipOpt=0;   
        checkChar=1;
    case 3
        fileName=varargin{1};
        T=varargin{2};
        skipOpt=varargin{3};
        checkChar=1;
    case 4
        fileName=varargin{1};
        T=varargin{2};
        skipOpt=varargin{3};
        checkChar=varargin{4};
end

% Force column if it isn't
if ~isvector(T)
    T=T(:);%Make column
end

% Convert non-character content to text
if checkChar==1
    indNotChar=find(~cellfun(@(x) ischar(x),T,'UniformOutput',1));
    if ~isempty(indNotChar)
        %warning('Non character entries detected');
        for q=indNotChar            
            T{q}=vec2strIntDouble(T{q},'%6.7e');
        end
    end
end

%%
% Write to file

fid=fopen(fileName,'w');
switch skipOpt
    case 0 %Process in one go
        fprintf(fid,'%s\r\n', T{:});
    case 1 %Looping to check if lines are empty
        for q=1:numel(T)
            l=T{q}; %Get current cell entry            
            if ~ischar(l) %If this isn't a char then attempt conversion
                l=sprintf('%u',l);                
            end
            
            %Write if not empty
            if ~isempty(deblank(l))
                fprintf(fid,'%s\r\n',l);
            end            
        end
end
fclose(fid);
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2023 Kevin Mattheus Moerman and the GIBBON contributors
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
