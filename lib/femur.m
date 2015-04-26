function [F,V]=femur

% [F,V]=femur
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% femur model. 
%
% This model consists of 5924 triangular faces and 2964 vertices.
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
fileName=fullfile(pathName,'femur.mat');
D=load(fileName);
F=D.F;
V=D.V;


