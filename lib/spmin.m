function [varargout]=spmin(varargin)

% function [minVal,minInd]=spmin(A,B,vecdim,nanflag,logicRelevant)
% ------------------------------------------------------------------------
%
%
% Kevin Moerman
% ------------------------------------------------------------------------

%%

A=varargin{1};

switch nargin
    case 5
        logicRelevant=varargin{5};
        varargin=varargin(1:4);
        if isempty(varargin{4})
            varargin=varargin(1:3);
        end
    otherwise
        logicRelevant=A~=0;
end

%%

maxOffset=max(A(logicRelevant));

A(logicRelevant)=A(logicRelevant)-maxOffset;
varargin{1}=A;
[minVal,minInd]=min(varargin{:});
minVal=minVal+maxOffset;

%% Collect output
varargout{1}=minVal;
varargout{2}=minInd;

