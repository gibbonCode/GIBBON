function struct2txt(varargin)

%% parse input

switch nargin
    case 2
        S=varargin{1};
        fileName=varargin{2};
        delimiter=[];
    case 3
        S=varargin{1};
        fileName=varargin{2};
        delimiter=varargin{3};
end

if isempty(delimiter)
    delimiter = ',';
end

%%

fileID = fopen(fileName,'w');

nameCell=fieldnames(S);
for q=1:1:numel(nameCell)
    entryNow=S.(nameCell{q});
    if isnumeric(entryNow)
        entryText=vec2strIntDouble(entryNow,'%6.7e');
    else
        entryText=entryNow;
    end
    textLine=[nameCell{q},delimiter,entryText];
    fprintf(fileID,[textLine,'\n']);
end
fclose(fileID);

