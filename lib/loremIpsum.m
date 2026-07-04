function [S]=loremIpsum(varargin)

% function [S]=loremIpsum(numWords,outputOpt)
% ------------------------------------------------------------------------
% This function creates the "lorem ipsum" test word set. It can be used as
% a simple text character set. The inputs include the number of words
% (numWords) requested and the output type (outputOpt) denoting wether the
% user requests a string or a cell. 
% 
% Change log: 
% 2023/09/01 KMM: Updated description/documentations
% ------------------------------------------------------------------------

%%

switch nargin
    case 0 
        numWords=69;
        outputOpt='string';
    case 1
        numWords=varargin{1};
        outputOpt='string';
    case 2
        numWords=varargin{1};
        outputOpt=varargin{2};
end

S='Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. ';
if numWords>69
    S=repmat(S,[1 ceil(numWords/69)]);
end
S=strsplit(S,' ')';
S=S(1:numWords);

switch outputOpt
    case 'string'
        S=char(S);
    case 'cell'
        
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
