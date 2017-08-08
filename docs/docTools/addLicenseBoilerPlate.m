clear; close all; clc;

%%

toolboxPath=fileparts(fileparts(fileparts(mfilename('fullpath'))));
docPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs');
libPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'lib');
docToolsPath=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'docs','docTools');

%%

licBoilerPlateFile=fullfile(toolboxPath,'licenseBoilerPlate.txt');
[T]=txtfile2cell(licBoilerPlateFile);

%%

replaceOn=1; 
fileExtension='.m';
targetStart='********** _license boilerplate_ **********';
