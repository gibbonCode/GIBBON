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
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2018  Kevin Mattheus Moerman
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
