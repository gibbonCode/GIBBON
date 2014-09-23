function [F,V]=subQuad(F,V,n)

%%
if size(V,2)==2
    V(:,3)=0;
end

for qIter=1:1:n
    
    numV=size(V,1);
    X=V(:,1); Y=V(:,2); Z=V(:,3); 

    %% DERIVE EDGE INDICES AND FACE-EDGE INDEX MATRIX
    
    %Format of column index in F
    EColumnInd=[(1:size(F,2)); (1:size(F,2))];
    EColumnInd=[EColumnInd(2:end) EColumnInd(1)];
    
    %Derive edges matrix
    E=F(:,EColumnInd)'; %Use index into F to create edges matrix
    E=reshape(E,2,numel(E)/2)';
        
    E=sort(E,2); %Sort edge order
    [E,~,ind2] = unique(E,'rows'); %Removing double edges, i.e. [1  4] = [4  1]      
    
    Fe=reshape(1:numel(F),size(F,2),size(F,1))';
    Fe=ind2(Fe);
    
    %% Calculate mid-face vertices
    
    XF=X(F); YF=Y(F); ZF=Z(F);
    V_midFace=[mean(XF,2) mean(YF,2) mean(ZF,2)];
    
    %% Calculate mid-edge vertices
    
    XE=X(E); YE=Y(E); ZE=Z(E);
    V_midEdge=[mean(XE,2) mean(YE,2) mean(ZE,2)];
    
    %% Create new faces matrix
    
    Vs=[V;V_midEdge;V_midFace];
    startIndMidEdge=numV;
    startIndMidFace=startIndMidEdge+size(V_midEdge,1);
    indMidFace=startIndMidFace+1:size(Vs,1);
    
    Fs=repmat(F,4,1);
    for q=1:1:4;
        startInd=1+(q-1)*size(F,1);
        endInd=startInd-1+size(F,1);
        
        if q==1
            ind4=4;
        else
            ind4=q-1;
        end
        fs=[F(:,q) Fe(:,q)+numV indMidFace(:) Fe(:,ind4)+numV];
        
        Fs(startInd:endInd,:)=fs;
    end   
   
    %% Overwrite faces and vertices
    V=Vs; F=Fs;   
    
end
