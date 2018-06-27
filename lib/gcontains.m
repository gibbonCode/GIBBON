function TF=gcontains(varargin)

try
    TF=contains(varargin{:}); %Added in R2018a
catch ME
    %Use alternative
    switch nargin
        case 2
            str=varargin{1};
            strPattern=varargin{2};
            ignoreValue=false;
        case {3,4}
            str=varargin{1};
            strPattern=varargin{2};
            % ignoreCase=varargin{3};
            ignoreValue=varargin{4};
        otherwise
            error('Wrong number of input arguments');
    end
    
    if ignoreValue
        if isa(str,'cell')
            for q=1:1:numel(str)
                str{q}=lower(str{q});
            end
        else
            str=lower(str);
        end
        
        if isa(strPattern,'cell')
            for q=1:1:numel(strPattern)
                strPattern{q}=lower(strPattern{q});
            end
        else
            strPattern=lower(strPattern);
        end
    end
    
    TF=false(size(str));
    if isa(strPattern,'cell')
        for q=1:1:numel(strPattern)
            TF=TF | ~cellfun(@isempty,strfind(str,strPattern{q}));
        end
    elseif isa(strPattern,'string')
        for q=1:1:numel(strPattern)
            TF=TF | ~cellfun(@isempty,strfind(str,strPattern{q}));
        end
    else
        TF=cellfun(@isempty,strfind(str,strPattern));
    end
    warning(ME.message)
end