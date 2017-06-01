function addHelpSearch

%% Make help files searchable

mfilePath = mfilename('fullpath');
gibbonPath=fileparts(fileparts(mfilePath));
gibbonHelpPath=fullfile(gibbonPath,'docs','html');
builddocsearchdb(gibbonHelpPath);
