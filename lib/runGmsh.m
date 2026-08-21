function runGmsh(geo_name, gmsh_path)

    if nargin < 2 || isempty(gmsh_path)

        if isunix
            execName = 'gmsh';
        else
            execName = 'gmsh.exe';
        end

        % Search in Additional Software (as installed by mpminstall),
        % and check the bare name, in case user has added it to path
        % gibbonSettings.set('GmshPath', PATH) will override both
        gmsh_path = {
            fullfile(gibbonSettings.gibbonPath, '..', 'Additional Software', 'Gmsh','**', execName);
            execName
        };
    end
    gmsh_path = gibbonSettings.findExec(gmsh_path, setting='GmshPath', fail='error');

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
