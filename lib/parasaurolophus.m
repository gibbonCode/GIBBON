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