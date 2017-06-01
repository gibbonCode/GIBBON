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
% 2017/06/01 Fixed bug in relation to operational system differences
%------------------------------------------------------------------------

%%
if ispc % semi-colon splits paths
    strPattern=[filesep,';'];    
else %colon splits paths for Linux and mac
    strPattern=':';        
end
pathNames=regexp(genpath(pathName),strPattern, 'split');
pathNames=pathNames(2:end-1)';

end