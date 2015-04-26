function [F,V]=utah_teapot

% [F,V]=utah_teapot
% ------------------------------------------------------------------------
%
% This function generates patch data (faces=F and vertices=V) for a
% utah_teapot model. 
%
% This model consists of 2464 triangular faces and 1232 vertices.
% 
% http://en.wikipedia.org/wiki/Utah_teapot
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
fileName=fullfile(pathName,'utah_teapot.mat');
D=load(fileName);
F=D.F;
V=D.V;


