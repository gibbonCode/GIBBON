function cleanUpTetGen(pathNameTempFiles)

% function cleanUpTetGen(pathNameTempFiles)
% ------------------------------------------------------------------------
% This function can clean a folder from tetgen output files. 
% 
% See also |cleanDir|
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------

extCell={'ele','node','face','edge','mtr','smesh','p2t'}; %Extensions of files to delete

%Remove the files with matching extensions
cleanDir(pathNameTempFiles,extCell);

