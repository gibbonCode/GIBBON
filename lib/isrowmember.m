function logicMember=isrowmember(varargin)

% function logicMember=isrowmember(F1,F2,orderOption)
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------

%% Parse input 

switch nargin    
    case 2
        F1=varargin{1};
        F2=varargin{2};
        orderOption=1;
    case 3
        F1=varargin{1};
        F2=varargin{2};
        orderOption=varargin{3};       
end

%%

%Sort if needed
if orderOption
    %Sort e.g. so that [4 3 2 1]=[1 2 3 4]
    F1=sort(F1,2);
    F2=sort(F2,2);
end

%Create virtual indicesf or each face
maxIndex=max([F1(:);F2(:)]);
sizVirt=maxIndex*ones(1,size(F1,2));
indVirt_F1=sub2indn(sizVirt,F1);
indVirt_F2=sub2indn(sizVirt,F2);

logicMember=ismember(indVirt_F1,indVirt_F2);

