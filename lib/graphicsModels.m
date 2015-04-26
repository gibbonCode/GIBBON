function [F,V]=graphicsModels(varargin)

% function [F,V]=graphicsModels(modelID)
% ------------------------------------------------------------------------
% 
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/04/25 Added to GIBBON
%------------------------------------------------------------------------

%% Parse input

switch nargin
    case 0
        modelID=1;
    case 1
        modelID=varargin{1};
end

%%
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');
switch modelID
    case 1
        fileName=fullfile(pathName,'stanford_bunny_closed.mat');
        D=load(fileName);
    case 2
        fileName=fullfile(pathName,'utah_teapot.mat');
        D=load(fileName);
    case 3
        fileName=fullfile(pathName,'cow.mat');
        D=load(fileName);
    case 4
        fileName=fullfile(pathName,'parasaurolophus.mat');
        D=load(fileName);
    case 5
        fileName=fullfile(pathName,'femur.mat');
        D=load(fileName);
    case 6
        fileName=fullfile(pathName,'hip_implant.mat');
        D=load(fileName);
end

F=D.F;
V=D.V;
