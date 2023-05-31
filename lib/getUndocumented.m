function getUndocumented(varargin)

% function getUndocumented(n)
% ------------------------------------------------------------------------
%
%
% ------------------------------------------------------------------------

%%
switch nargin
    case 0
        n=1;
    case 1
        n=varargin{1};
end

%%
testFolder=fullfile(fileparts(fileparts(mfilename('fullpath'))),'docs');
[testFileList]=getTestFiles('HELP'); 

c=1; 
for q=1:1:numel(testFileList)

    mFileNow=fullfile(testFolder,testFileList{q});

    % fid = fopen(mFileNow);
    % lineNow = fgetl(fid) 

    fileLines = readlines(mFileNow);

    if any(contains(fileLines,'UNDOCUMENTED'))
        if c==n
            open(mFileNow)
            return
        end
        c=c+1;
    end    
end

