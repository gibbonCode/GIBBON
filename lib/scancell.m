function [IND_found]=scancell(T,targets,found_count)

line_count=1;
IND_found=cell(numel(targets),1);

if isempty(found_count)
    while 1
        for i=1:1:numel(targets)
            if ~isempty(strfind(T{line_count},targets{i}))
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
            if ~isempty(strfind(T{line_count},targets{i}))
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
