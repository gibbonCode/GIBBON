function [S]=txt2struct(varargin)

%% parse input

switch nargin
    case 1
        fileName=varargin{1};
        delimiter=[];
    case 2
        fileName=varargin{1};
        delimiter=varargin{2};
end

if isempty(delimiter)
    delimiter = ',';
end
%%
if exist(fileName,'file')
    %Import data into cell
    
    formatSpec = '%s%[^\n\r]';
    fileID = fopen(fileName,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    
    % Convert cell to structure        
    nameCell=dataArray{1};
    entryCell=dataArray{2};    
    [S]=cellPair2struct(nameCell,entryCell,1);
end