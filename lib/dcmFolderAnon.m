function dcmFolderAnon(pathName,varargin)

% function dcmFolderAnon(pathName,varargin)
% ------------------------------------------------------------------------
%
% This function anonomizes the dicom files in the folder pathName using the
% dicomanon function. The first input is the folder path name, the
% following inputs are possible inputs of the dicomanon function. 
%
% See also: dicomanon
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2018/02/09: Created
%------------------------------------------------------------------------

%%

files = dir(fullfile(pathName,'*.dcm'));
files={files(1:end).name};
files=sort(files(:));

hw = waitbar(0,'Anonomizing DICOM info...');  
for q=1:1:numel(files)   
    fileName=fullfile(pathName,files{q});
    dicomanon(fileName,fileName,varargin{:});
    waitbar(q/numel(files),hw,['Anonomizing DICOM info...',num2str(round(100.*q/numel(files))),'%']);
end
close(hw)