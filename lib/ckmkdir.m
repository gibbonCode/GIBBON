function [varargout]=ckmkdir(dirName)

%%

% Check for existence of directory
existFlag=exist(dirName,'dir');
if ~existFlag %Check if folder exists   
    mkdir(dirName); %Make the directory if it does not exist
end

%Store output
if nargout==1
    varargout{1}=existFlag;
end