function [F,V]=elephant

% function [F,V]=elephant
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for an
% elephant model.
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/10/22
%------------------------------------------------------------------------
%%

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
offPath=fullfile(toolboxPath,'data','OFF');
fileName=fullfile(offPath,'elephant-50kv.off');
[F,V] = import_off(fileName);


