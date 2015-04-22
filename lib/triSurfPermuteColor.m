function [C,logicConverged]=triSurfPermuteColor(varargin)


%% Parse input

switch nargin 
    case 2
        F=varargin{1};
        V=varargin{2};        
        nc=4;
        maxIter=1000;
    case 3
        F=varargin{1};
        V=varargin{2};        
        nc=varargin{3};
        maxIter=1000;
    case 4
        F=varargin{1};
        V=varargin{2};
        nc=varargin{3};
        maxIter=varargin{4};
    otherwise
        error('Wrong number of input arguments');
end

C=randi(nc,size(F,1),1);

%%

TR = triangulation(F,V);
N = neighbors(TR);

%%
numFaces=size(F,1);

%% Random iterations

q=0;
logicConverged=0;
numPrev=numFaces;
while q<maxIter && logicConverged==0
    
    q=q+1;
    
    CN=C(N); %Colors for neighbours
    
    logicDoubleColor=CN==C(:,ones(size(CN,2),1));    
    logicRep=any(logicDoubleColor,2);
    
    C(logicRep)=randi(nc,nnz(logicRep),1);
    
    %     set(hp,'CData',C);
    %     drawnow;
    numNow=nnz(logicRep);
    p=abs(numPrev-numNow)./numFaces;
    
    numPrev=nnz(logicRep);
    
    logicConverged=nnz(logicRep)==0;
    
    if q==maxIter
        break
    end 
end

logicConverged=nnz(logicRep)==0;

%% Structured iterations

if ~logicConverged
    while logicConverged==0
        
        CN=C(N); %Colors for neighbours
        
        logicDoubleColor=CN==C(:,ones(size(CN,2),1));
        logicRep=any(logicDoubleColor,2);
        
        Ir=find(logicRep,1);
        
        c=C(Ir);
        cn=CN(Ir,:);
        
        cFound=0;
        s=1;
        while cFound
            if ~any(cn-s==0);
                C(Ir)=s;
                cFound=1;
            elseif s==nc
                cFound=1;
            end
            s=s+1;
        end
        %         set(hp,'CData',C);
        %         drawnow;
        
        if nnz(logicRep)==0
            logicConverged=1;
        else
            q=q+1;
        end
        
    end
end
logicConverged=nnz(logicRep)==0;

