function [F,V]=parasaurolophus

% [F,V]=parasaurolophus
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% triangulated parasaurolophus dinosaur model. 
%
% This adjusted model consists of 892 triangular faces and 448 vertices. 
% 
% The model was constructed based on the model given here: 
% https://www.rocq.inria.fr/gamma/gamma/download/affichage.php?dir=DINOSAUR/&name=Parasaurolophus
% 
%
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2015/04/25 Updated for GIBBON
%------------------------------------------------------------------------
%%

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');
fileName=fullfile(pathName,'parasaurolophus.mat');
D=load(fileName);
F=D.F;
V=D.V;
 
%% <-- GIBBON footer text --> 
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
