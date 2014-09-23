function runTetView(modelName)

%% SETTING TETGEN PATHNAMES
% cdNow=pwd;
pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen');
runNameTetView=fullfile(pathNameTetView,'tetview-win.exe');

modelName=regexprep(modelName,'\','/');

%% RUN TETVIEW
% cd(pathNameTetView);

runString=[runNameTetView,' ',modelName,' & '];
system(runString);

% cd(cdNow);

