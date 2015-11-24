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

% filePath=mfilename('fullpath');
% toolboxPath=fileparts(fileparts(filePath));
% offPath=fullfile(toolboxPath,'data','OFF');
% fileName=fullfile(offPath,'elephant-50kv.off');
% [F,V] = import_off(fileName);
% 
% [F,V,~,~]=triSurfRemoveThreeConnect(F,V,[]);
% 
% [R,~]=euler2DCM([1/3*pi 0 0]);
% V=(R*V')';
% [R,~]=euler2DCM([0 1/36*pi 0]);
% V=(R*V')';

[F,V]=graphicsModels(7);

