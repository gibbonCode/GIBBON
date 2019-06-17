function TF=gcontains(varargin)

try    
    TF=contains(varargin{:}); %Added in R2018a
catch ME
    %Use alternative
    switch nargin
        case 2
            str=varargin{1};
            strPattern=varargin{2};
            ignoreValue=false;
        case {3,4}
            str=varargin{1};
            strPattern=varargin{2};
            % ignoreCase=varargin{3};
            ignoreValue=varargin{4};
        otherwise
            error('Wrong number of input arguments');
    end
    
    if ignoreValue
        if isa(str,'cell')
            for q=1:1:numel(str)
                str{q}=lower(str{q});
            end
        else
            str=lower(str);
        end
        
        if isa(strPattern,'cell')
            for q=1:1:numel(strPattern)
                strPattern{q}=lower(strPattern{q});
            end
        else
            strPattern=lower(strPattern);
        end
    end
    
    TF=false(size(str));
    if isa(strPattern,'cell')
        for q=1:1:numel(strPattern)
            TF=TF | ~cellfun(@isempty,strfind(str,strPattern{q}));
        end
    elseif isa(strPattern,'string')
        for q=1:1:numel(strPattern)
            TF=TF | ~cellfun(@isempty,strfind(str,strPattern(q)));
        end
    elseif isa(strPattern,'char')        
        TF=strfind(str,strPattern); 
        if isa(TF,'cell')
            TF=~cellfun(@isempty,TF);
        else
            TF=~isempty(TF);
        end        
    else
        TF=~cellfun(@isempty,strfind(str,strPattern));
    end
    % warning(ME.message)
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
