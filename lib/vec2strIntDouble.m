function [t]=vec2strIntDouble(varargin)

%%
switch nargin
    case 1
        n=varargin{1};
        formatLong='%6.7e';
    case 2
        n=varargin{1};
        formatLong=varargin{2};
end

%%

if isnumeric(n) %If it is numeric
    n=double(n);
    if isrounded(n) %If it looks like an integer
        t_form=repmat('%d, ',1,size(n,2)); 
    else %Not an integer
        t_form=repmat([formatLong,', '],1,size(n,2)); 
    end
    t_form=t_form(1:end-2); %Take away last space and comma
    
    %Convert to string
    t=sprintf(t_form,n);
elseif char(n)
    t=n;
else
    error('Input should be numeric')
end