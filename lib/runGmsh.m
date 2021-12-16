
function runGmsh(geo_name,gmsh_path)
    if isempty(gmsh_path)
        gmsh_path='C:\Users\icberg\OneDrive - UAB - The University of Alabama at Birmingham\Documents\FEM Model Heart Chip';
    end
    
    oldFolder = cd(gmsh_path);
    system(['gmsh.exe ' geo_name ' -0']);
    cd(oldFolder)
end    