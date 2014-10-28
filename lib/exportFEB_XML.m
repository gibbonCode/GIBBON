function exportFEB_XML(varargin)

%% Parse input
switch nargin
    case 2
        save_name=varargin{1};
        XDOC=varargin{2};
        topCommentLine=['Created using GIBBON, ',datestr(now)]; %Add GIBBON comment at start
    case 3
        save_name=varargin{1};
        XDOC=varargin{2};
        topCommentLine=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

%%
%Write to text file using xmlwrite (adds extra lines for unknown reasons)
xmlwrite(save_name,XDOC);

%Import back into cell array
[T]=txtfile2cell(save_name);

if ~isempty(topCommentLine)
    TT=cell(1,1);
    TT(1,1)=T(1,1);
    TT(2,1)={'<!-- '};
    TT(3,1)={topCommentLine};
    TT(4,1)={'-->'};
    TT(end+1:end+size(T,1)-1)=T(2:end);
else
    TT=T; 
end

%Now save to txt file while skipping "empty lines"
cell2txtfile(save_name,TT,1);


