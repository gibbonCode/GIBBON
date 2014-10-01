function [varargout]=maxnumel(a)

%Max variable size available
[numFreeBytes]=freeMemory;

%Size of input variable
sizA=size(a);
whos_struct=whos('a');
numBytesType=whos_struct.bytes;

%Max number of input type available
n=floor(numFreeBytes./numBytesType);

varargout{1}=n;
switch nargout    
    case 2
        varargout{2}=numBytesType;
end
    
end