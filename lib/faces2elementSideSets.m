function [varargout]=faces2elementSideSets(varargin)
% function [elementId,elementSides,elemSets,elemSetSides]=faces2elementSideSets(f,E)
% ------------------------------------------------------------------------
%
% Kevin Mattheus Moerman
% 2020/05/30 Created
% ------------------------------------------------------------------------

%%

switch nargin
    case 2
        f=varargin{1};
        E=varargin{2};
        reorderOpt=0;
    case 3
        f=varargin{1};
        E=varargin{2};
        reorderOpt=varargin{3};
end
%%
logicTouch=sum(ismember(E,f),2)==size(f,2);
ind=(1:1:size(E,1))';
[Fs,C,CF]=element2patch(E(logicTouch,:),ind(logicTouch));
logicMember=isrowmember(Fs,f,reorderOpt);
elementId=C(logicMember);
elementSides=CF(logicMember);
varargout{1}=elementId;
varargout{2}=elementSides;
if nargout>2
    setSideUni=unique(elementSides);
    numSets=numel(setSideUni);
    elemSets=cell(1,numSets);
    elemSetSides=zeros(1,numSets);
    for q=1:1:numSets
        sideType=setSideUni(q);
        logicNow=elementSides==sideType;
        indElem=unique(elementId(logicNow));
        elemSets{q}=indElem(:)';
        elemSetSides(q)=sideType;
    end
    varargout{3}=elemSets;
    varargout{4}=elemSetSides;
end

end