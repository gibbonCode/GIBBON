function runGmsh(geo_name, gmsh_path)

    if nargin < 2 || isempty(gmsh_path)

        % as downloaded when installing GIBBON with mpminstall
        gmsh_path = fullfile(gibbonSettings.gibbonPath, '..', 'Additional Software', 'Gmsh');

        if isunix
            gmsh_path = fullfile(gmsh_path,'**/bin/gmsh');
        else
            gmsh_path = fullfile(gmsh_path, '**', 'gmsh.exe');
        end
        d = dir(gmsh_path);
        if isempty(d)
            error('GMSH executable not found');
        end
        gmsh_path = fullfile(d.folder, d.name);
    end

    cmd = sprintf('"%s" "%s" -0', gmsh_path, geo_name);
    status = system(cmd);
    assert(status == 0, 'GMSH execution failed, used command:\n%s', cmd);
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
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
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
