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
