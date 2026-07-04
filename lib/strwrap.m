function str=strwrap(varargin)

% function str=strwrap(str,n,pattern,replaceOpt)
%-------------------------------------------------------------------------
% This function wraps the input string 
% Example, if the input is: 
%   str='1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18';
%   n=8;
%   pattern=', '
% Then the output is: 
% str=    '1, 2, 3, 4, 5, 6, 7, 8, 
%          9, 10, 11, 12, 13, 14, 15, 
%          16, 17, 18'
% I.e. a new line is created every n pattern seperated entries. 
%
% Change log: 
% 2023/03/10: The default behaviour is now set to add a new line character,
% rather than to replace the pattern by a new line character. The option to
% switch between these two has been added. 
% 2023/04/20: Improved speed of replacement of nth occurance by newline
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        str=varargin{1};
        n=[];
        pattern=[];
        replaceOpt=false;
    case 2
        str=varargin{1};
        n=varargin{2};
        pattern=[];
        replaceOpt=false;
    case 3
        str=varargin{1};
        n=varargin{2};
        pattern=varargin{3};
        replaceOpt=false;
    case 4            
        str=varargin{1};
        n=varargin{2};
        pattern=varargin{3};
        replaceOpt=varargin{4};
end

if isempty(n)
    n=1;
end

if isempty(pattern)
    pattern=', ';
end

%%

% try %Use count as introduced in R2016b
%     N = count(str,pattern); 
% catch %Use old approach
%     N = numel(strfind(str,pattern));
% end
% 
% rangeSteps=n:n:N;
% 
% %Replace certain occurances of the pattern with a new line symbol
% for q=1:1:numel(rangeSteps)    
%         str=regexprep(str,pattern,'\n',rangeSteps(q)-(q-1));    
% end

ind=strfind(str,pattern);
str(ind(n:n:end))=newline;

%Put the pattern back at these locations if needed
if replaceOpt==false
    str=regexprep(str,'\n',[pattern,'\n']);    
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
