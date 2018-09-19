function str=strwrap(varargin)

% function str=strwrap(str,n,pattern)
%-------------------------------------------------------------------------
%
%-------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        str=varargin{1};
        n=[];
        pattern=[];
    case 2
        str=varargin{1};
        n=varargin{2};
        pattern=[];
    case 3
        str=varargin{1};
        n=varargin{2};
        pattern=varargin{3};
end

if isempty(n)
    n=1;
end

if isempty(pattern)
    pattern=', ';
end

%%
N = count(str,pattern);
rangeSteps=n:n:N;
for q=1:1:numel(rangeSteps)
    str=regexprep(str,pattern,'\n',rangeSteps(q)-(q-1));
end
