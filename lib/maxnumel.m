function [varargout]=maxnumel(a)

%Max variable size available
mem_struct=memory;
num_bytes=mem_struct.MaxPossibleArrayBytes;

%Size of input variable
whos_struct=whos('a');
num_bytes_type=whos_struct.bytes;

%Max number of input type available
n=floor(num_bytes./num_bytes_type);

varargout{1}=n;
switch nargout    
    case 2
        varargout{2}=num_bytes_type;
end
    
end