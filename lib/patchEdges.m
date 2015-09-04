function E=patchEdges(varargin)

% function E=patchEdges(varargin)
% -----------------------------------------------------------------------
% E=patchEdges(F,uniOpt)
% Uses the input faces array F to compute an edge array E. If uniOpt==1
% then the output array contain unique edges (irrespective of node order
% such that e.g. [4 1] and [1 4] are seen as the same edge). If uniOpt~=1
% then the double edges are maintained. If only one input is provided it is
% assumed to represent F and the default, uniOpt=0, is used. 
%
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2014/03/17
%------------------------------------------------------------------------

%% PARSE INPUT
F=varargin{1};
switch nargin
    case 1
        uniOpt=0;
    case 2
        uniOpt=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end
%% DERIVE NON-UNIQUE EDGES MATRIX
E1=F';
E2=F(:,[2:end 1])';
E=[E1(:) E2(:)];

%% REMOVE DOUBLE ENTRIES IF DESIRED

if uniOpt==1
    Es=sort(E,2); %Sorted so [1 4] and [4 1] are seen as the same edge
    [~,indUni,~]=unique(Es,'rows'); %Get indices for unique edges
    E=E(indUni,:);
end

