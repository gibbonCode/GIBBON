function S=inputdlgStruct(varargin)

%% Parse input
switch nargin
    case 1        
        nameCell=varargin{1};
        defaultOptions=[];
        dialogTitle=[];
    case 2        
        nameCell=varargin{1};
        defaultOptions=varargin{2};
        dialogTitle=[];
    case 3
        nameCell=varargin{1};
        defaultOptions=varargin{2};
        dialogTitle=varargin{3};
end

if isempty(dialogTitle)
   dialogTitle='Input dialog'; 
end

if isempty(defaultOptions)
   defaultOptions=repmat({''},size(nameCell));
end
  
%%
s=25+max([cellfun(@numel,nameCell) cellfun(@numel,defaultOptions)]); %Set sizes
entryCell = inputdlg(nameCell,dialogTitle,[1 s],defaultOptions); %Open dialog box
[S]=cellPair2struct(nameCell,entryCell,1); %Convert output to structure