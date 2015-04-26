function [F,V]=hip_implant

% [F,V]=hip_implant
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% hip_implant model. 
%
% This model consists of 20958 triangular faces and 10483 vertices.
% 
%
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2015/04/26
%------------------------------------------------------------------------
%%

defaultFolder = fileparts(fileparts(mfilename('fullpath')));
pathName=fullfile(defaultFolder,'data','libSurf');
fileName=fullfile(pathName,'hip_implant.mat');
D=load(fileName);
F=D.F;
V=D.V;


