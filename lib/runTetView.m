function [varargout]=runTetView(modelName)

%% SETTING TETGEN PATHNAMES

compString=computer; 
switch compString
    case 'PCWIN' %Windows 32-bit
        error('PCWIN 32-bit is not supported. Compile tetGen from the source and alter the code here');
    case 'PCWIN64' %Windows 64-bit
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','win64');
        runNameTetView=fullfile(pathNameTetView,'tetview-win.exe');
    case 'GLNXA64'        
        pathNameTetView=fullfile(fileparts(fileparts(mfilename('fullpath'))),'lib_ext','tetGen','lin64');
        runNameTetView=fullfile(pathNameTetView,'tetview-linux');
    case 'MACI64'        
        error('MACI64 is not supported yet. Get TetView online and alter the code here');
    otherwise
        error('Your platform does not seem to be supported. Code your own solution or contact support.')
end

modelName=regexprep(modelName,'\','/');

%% RUN TETVIEW

runString=['"',runNameTetView,'" "',modelName,'" & '];
[runStatus,runCmdHist]=system(runString);

varargout{1}=runStatus;
varargout{2}=runCmdHist;

