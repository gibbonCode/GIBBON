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
    case 7 
        fileName=fullfile(pathName,'elephant.mat');
        D=load(fileName);
end

F=D.F;
V=D.V;
 
%% 
% ********** _license boilerplate_ **********
% 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
