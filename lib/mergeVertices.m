function [varargout]=mergeVertices(varargin)

% function [F,V,ind1,ind2]=mergeVertices(F,V,nKeep)


%% Parse input

switch nargin    
    case 2
        F=varargin{1};
        V=varargin{2};
        nKeep=5;
    case 3
        F=varargin{1};
        V=varargin{2};
        nKeep=varargin{3};
    otherwise
        error('Wrong number of input arguments');
end

%% Merge nodes

[~,ind1,ind2]=unique(pround(V,nKeep),'rows');
V=V(ind1,:);
F=ind2(F);

%% Collect output

varargout{1}=F;
varargout{2}=V;
varargout{3}=ind1;
varargout{4}=ind2;

%%