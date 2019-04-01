function [varargout]=spmax(varargin)

% function [maxVal,maxInd]=spmax(varargin)
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

minOffset=min(A(logicRelevant));

A(logicRelevant)=A(logicRelevant)-minOffset;
varargin{1}=A;
[maxVal,maxInd]=max(varargin{:});
maxVal=maxVal+minOffset;

%% Collect output
varargout{1}=maxVal;
varargout{2}=maxInd;

