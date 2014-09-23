function [pathNames]=getSubPaths(pathName)

% function [pathNames]=getSubPaths(pathName)
% ------------------------------------------------------------------------
%
%This function (based on the GENPATH command) creates the output pathNames
%which contains all folders and sub-folders within the folder specified by
%the input pathName.
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 18/04/2013
%------------------------------------------------------------------------

pathNames=regexp(genpath(pathName), '\;', 'split');
pathNames=pathNames(2:end-1)';

end