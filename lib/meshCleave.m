function [varargout]=meshCleave(varargin)

% function [logicAt,logicAbove,logicBelow]=meshCleave(E,V,P,n,inclusiveSwitch)
% ------------------------------------------------------------------------
%
% 
% Kevin Mattheus Moerman
% 2020/05/18: Created
%
% ------------------------------------------------------------------------

%% Parse input
switch nargin
    case 2
        E=varargin{1};
        V=varargin{2};        
        P=mean(V,1);
        n=[0 0 1]; 
        inclusiveSwitch=[0 0];
    case 3
        E=varargin{1};
        V=varargin{2};        
        P=varargin{3};
        n=[0 0 1];        
        inclusiveSwitch=[0 0];
    case 4
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};        
        inclusiveSwitch=[0 0];
    case 5
        E=varargin{1};
        V=varargin{2};
        P=varargin{3};
        n=varargin{4};
        inclusiveSwitch=varargin{5};
end

%Check P
if isempty(P)
    P=mean(V,1); %Use default as mean of coordinates
end

%Check normal direction
if isempty(n)
    n=[0 0 1]; %Default z-direction
end
n=vecnormalize(n); %Normalize n vector

%% Construct rotation matrix
nz=[0 0 1]; %The z-vector
R=vecPair2Rot(n,nz);

%%

%Rotate coordinates
Vr=V-P(ones(size(V,1),1),:);
Vr=Vr*R';
Z=Vr(:,3);

%%

if inclusiveSwitch(1)==1
    logicVerticesBelow=Z<=0; %Logic for nodes below slice
else
    logicVerticesBelow=Z<0; %Logic for nodes below slice
end

if inclusiveSwitch(2)==1
    logicVerticesAbove=Z>=0; %Logic for nodes above slice
else
    logicVerticesAbove=Z>0; %Logic for nodes above slice
end

logicAt=any(logicVerticesBelow(E),2) & any(logicVerticesAbove(E),2); %Logic for slice elements
logicBelow=all(logicVerticesBelow(E),2); %Logic for elements with nodes all below
logicAbove=all(logicVerticesAbove(E),2); %Logic for elements with nodes all above

%%

varargout{1}=logicAt;
varargout{2}=logicAbove;
varargout{3}=logicBelow;

end