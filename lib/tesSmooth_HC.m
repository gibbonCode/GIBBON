function [P]=tesSmooth_HC(TES,V,IND_V,cPar)

%% CONTROL PARAMETERS

if isfield(cPar,'Alpha')
    alp=cPar.Alpha;
else
    alp=0.1; %DEFAULT
end

if isfield(cPar,'Beta')
    bet=cPar.Beta;
else
    bet=0.5; %DEFAULT
end

if isfield(cPar,'n')
    nMax=cPar.n;
else
    nMax=1; %DEFAULT
end

if isfield(cPar,'RigidConstraints')
    indRigid=cPar.RigidConstraints;
else
    indRigid=[]; %DEFAULT
end

if isfield(cPar,'Tolerance')
    SSQD_Tol=cPar.Tolerance;
else
    SSQD_Tol=[]; %DEFAULT
end

%%

if ~isempty(SSQD_Tol)
    SSQD_old=[];
    SSQD_ratio=0;
end

if isempty(IND_V)
    [~,IND_V]=patchIND(TES,V,2);
end
logicValid=IND_V>0;

%%

nDims=size(V,2); %Number of dimensions

VP=NaN(size(IND_V,1),size(IND_V,2),nDims);

%Initializing coordinate sets
O=V; %Original point set
P=V; %"New" point set
Q=V; %"Previous" point set
B=V; %
C=V; %

for qIter=1:nMax;   
        
    %% SMOOTHENING
    
    %Loop for all dimensions
    for qDim=1:1:nDims
        
        %Simple Laplacian operation
        Xq=VP(:,:,qDim);
        Xq(logicValid)=Q(IND_V(logicValid),qDim);
        Xq=nanmean(Xq,2);           
        P(:,qDim)=Xq;   
        
        %
        B(:,qDim)=P(:,qDim)-((alp*O(:,qDim))+(((1-alp)*Q(:,qDim))));
        
        %
        Xb=VP(:,:,qDim);
        Xb(logicValid)=B(IND_V(logicValid),qDim);
        Xb=nanmean(Xb,2);
        P(:,qDim)= P(:,qDim)-((bet*B(:,qDim))+(((1-bet)*Xb)));
        
    end
        
    %%
        
    %Put back constrained points
    if ~isempty(indRigid)
       P(indRigid,:)=V(indRigid,:);
    end
    
    if ~isempty(SSQD_Tol)
        %Compute sum of squared differences with respect to previous iteration
        SSQD_new=nansum((P(:)-Q(:)).^2);
        if ~isempty(SSQD_old)
            SSQD_ratio=SSQD_new./SSQD_old;            
        end
        
        %Store current metrics        
        SSQD_old=SSQD_new;        
        if abs(1-SSQD_ratio)<=SSQD_Tol
            disp('Tolerance reached!');
           break %STOP SMOOTHING LOOP IF TOLERANCE IS REACHED 
        end
    end
    Q=P; %Store current metric
    
end

