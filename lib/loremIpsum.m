function S=loremIpsum(varargin)

switch nargin
    case 0 
        numWords=69;
        outputOpt='string';
    case 1
        numWords=varargin{1};
        outputOpt='string';
    case 2
        numWords=varargin{1};
        outputOpt=varargin{2};
end

S='Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum. ';
if numWords>69
    S=repmat(S,[1 ceil(numWords/69)]);
end
S=strsplit(S,' ')';
S=S(1:numWords);

switch outputOpt
    case 'string'
        S=char(S);
    case 'cell'
        
end