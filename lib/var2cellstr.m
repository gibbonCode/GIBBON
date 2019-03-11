function T=var2cellstr(varargin)

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

%Fix option structure, complete and remove empty values
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,1); %Complement provided with default if missing or empty

%Get parameters
formatDouble=optionStruct.formatDouble;
formatInteger=optionStruct.formatInteger;
dlmChar=optionStruct.dlmChar;

if isnumeric(A)
    n=size(A,2);
    A=A'; %transpose
    A=double(A); %Convert to double
    
    %Create string convertion format
    if all(isrounded(A)) %If it looks like an integer
        t_form=repmat([formatInteger,dlmChar,' '],1,n);
    else %Not an integer
        t_form=repmat([formatDouble,dlmChar,' '],1,n);
    end
    t_form=t_form(1:end-1-numel(dlmChar)); %Take away last space and comma
    t_form=[t_form,'\n']; %Append end of line character
    
    %Convert to string
    t=sprintf(t_form,A);
    t=t(1:end-1); %take away last end of line symbol
    T=splitlines(t);
else
    switch class(A)
        case 'cell'
            logicNumeric=cellfun(@(x) isnumeric(x),A);
            T=repmat({''},size(A));
            if any(logicNumeric)
                A_numeric=A(logicNumeric);
                numEntries=cellfun(@(x) numel(x),A_numeric);
                T_numeric=repmat({''},size(A_numeric));
                if any(numEntries>1)
                    for q=1:1:numel(A_numeric)
                        a=A_numeric{q};
                        T_numeric{q}=mat2strIntDouble(a);
                    end
                else
                    logicRounded=cellfun(@(x) isrounded(x),A_numeric);
                    
                    if any(logicRounded)
                        A_numeric_rounded=A_numeric(logicRounded);
                        t = sprintf('%d \n',A_numeric_rounded{:});
                        t=t(1:end-1); %take away last end of line symbol
                        T_numeric_rounded=splitlines(t);
                        
                        A_numeric_notrounded=A_numeric(~logicRounded);
                        t = sprintf('%6.7e \n',A_numeric_notrounded{:});
                        t=t(1:end-1); %take away last end of line symbol
                        T_numeric_notrounded=splitlines(t);
                        
                        T_numeric(logicRounded)=T_numeric_rounded;
                        T_numeric(~logicRounded)=T_numeric_notrounded;
                    else
                        t = sprintf('%6.7e \n',A_numeric{:});
                        t=t(1:end-1); %take away last end of line symbol
                        T_numeric=splitlines(t);
                    end
                end
                T(logicNumeric)=T_numeric;
                T(~logicNumeric)=A(~logicNumeric);
            else
                T=A;
            end
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
