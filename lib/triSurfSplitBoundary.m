function [varargout]=triSurfSplitBoundary(varargin)

% function [F,V]=triSurfSplitBoundary(F,V,E,n,CF,CV)
% ------------------------------------------------------------------------
%
% 
% To do: Use for loop and recursion instead
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 3
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=size(E,1)+1;
        CF=[];
        CV=[];
    case 4
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=[];
        CV=[];
    case 5
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=varargin{5};
        CV=[];
    case 6
        F=varargin{1};
        V=varargin{2};
        E=varargin{3};
        n=varargin{4};
        CF=varargin{5};
        CV=varargin{6};
end

%%
if isempty(E)
    error('Edge set is empty');
end

if size(E,1)==n
    %Do nothing: Input = output
elseif size(E,1)>n
    error('Current number of edges exceeds number of desired edges');
else
    while size(E,1)~=n
        [~,indMax]=max(edgeLengths(E,V));
        numNodesPre=size(V,1); %Number of nodes before split
        [F,V,~,CF,CV]=triEdgeSplit(F,V,E(indMax,:),CF,CV);
        numNodesPost=size(V,1); %Number of nodes after split
        Eb=patchBoundary(F,V); %All boundary edges
        logicKeep=all(ismember(Eb,E) |...
            ismember(Eb,(numNodesPre+1):numNodesPost),2);
        E=Eb(logicKeep,:);
    end
end

%% Collect output
varargout{1}=F;
varargout{2}=V;
varargout{3}=E;
varargout{4}=CF;
varargout{5}=CV;

