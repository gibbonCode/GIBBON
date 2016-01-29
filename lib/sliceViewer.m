function hf=sliceViewer(varargin)


%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
        viewerType=1;
    case 2
        M=varargin{1};
        v=varargin{2};
        viewerType=1;
    case 3
        M=varargin{1};
        v=varargin{2};
        viewerType=varargin{3};
end
        
%% Start viewer

switch viewerType
    case 1
        hf=sv2(M,v);
    case 2
        hf=sv3(M,v);
end

