function [t]=mat2strIntDouble(varargin)

% function [t]=mat2strIntDouble(A,optionStruct)
%-------------------------------------------------------------------------
%
% Change log; 
% 2018/09/06 Altered behavior so column vectors are not forced to row
% vector text output. 
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        A=varargin{1};
        optionStruct=[];
    case 2
        A=varargin{1};
        optionStruct=varargin{2};
end

%Default structure
defaultOptionStruct.formatDouble='%6.7e';
defaultOptionStruct.formatInteger='%d';
defaultOptionStruct.dlmChar=',';
% defaultOptionStruct.forceRow=1;

%Fix option structure, complete and remove empty values
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%Get parameters
formatDouble=optionStruct.formatDouble;
formatInteger=optionStruct.formatInteger;
dlmChar=optionStruct.dlmChar;

%%
  
% Create character string
if isnumeric(A) %If it is numeric
    
    % Alter behaviour based on vector/matrix input
    n=size(A,2);
    A=A'; %transpose
    
    A=double(A); %Convert to double
    if all(isrounded(A)) %If it looks like an integer
        t_form=repmat([formatInteger,dlmChar,' '],1,n);        
    else %Not an integer
        t_form=repmat([formatDouble,dlmChar,' '],1,n); 
    end
    t_form=t_form(1:end-1-numel(dlmChar)); %Take away last space and comma

    %Append end of line character
    t_form=[t_form,'\n'];
    
    %Convert to string
    t=sprintf(t_form,A);
    
    %Take away last end of line character
    t=t(1:end-1);

elseif ischar(A)
    t=A;
elseif iscell(A)
    if ischar([A{:}])
        t = sprintf('%s \n',A{:}); 
        t=t(1:end-1);
    else
        A=reshape([A{:}],size(A));
        [t]=mat2strIntDouble(A,optionStruct);
    end
else
    error('Input should be numeric')
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
