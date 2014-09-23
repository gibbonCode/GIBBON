function varargout=minDist(varargin)

% function [D1,minIND]=minDist(V1,V2,maxVarSize)

%% Parse input
if nargin<2
    error('Insufficient input arguments');
end

V1=varargin{1};
V2=varargin{2};

if nargin>=3
    maxVarSize=varargin{3};
else
    %Max variable size available
    mem_struct=memory;
    num_bytes=mem_struct.MaxPossibleArrayBytes;    
    maxVarSize=num_bytes/2;    
end

%Derive class dependent variable size
[~,b1]=maxnumel(V1(1));
[~,b2]=maxnumel(V2(1));
b=max([b1 b2]);
numelVar=numel(V1)*numel(V2);
varSize=numelVar*b;

numSteps=ceil(varSize/maxVarSize);
indSteps=round(linspace(0,size(V1,1),numSteps));
indSteps=sort(unique(indSteps));
numSteps=numel(indSteps);

if numSteps>1 %In steps
    D1=zeros(size(V1,1),1);
    minIND=zeros(size(V1,1),1);
    for q=1:1:numSteps-1
        v1=V1(indSteps(q)+1:indSteps(q+1),:);
        try 
            d=dist(v1,V2'); %dist from Neural network toolbox
        catch
            d=distND(v1,V2); %GIBBON's dist function
        end
        [min_d,min_ind]=min(d,[],2);
        D1(indSteps(q)+1:indSteps(q+1))=min_d;
        minIND(indSteps(q)+1:indSteps(q+1))=min_ind;        
    end
else %In one go
    try
        D=dist(V1,V2'); %dist from Neural network toolbox
    catch
        D=distND(V1,V2); %GIBBON's dist function
    end
    [D1,minIND]=min(D,[],2);         
    D1=D1(:);
    minIND=minIND(:);
end

switch nargout
    case 1
        varargout{1}=D1;
    case 2
        varargout{1}=D1;
        varargout{2}=minIND;
    otherwise
        error('wrong number of output arguments');
end
