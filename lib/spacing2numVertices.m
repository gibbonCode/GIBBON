function nP=spacing2numVertices(F,V,pointSpacing)

A=sum(patch_area(F,V)); %Total area
At=(pointSpacing.^2*sqrt(3))/4; %Theoretical area of equilateral triangle
NF=(A/At);

E=patchEdges(F,1); %Edges
nE=size(E,1); %Number of edges
nF=size(F,1); %Number of faces
nV=size(V,1); %Number of vertices

%Compute Euler characteristic 
X=size(V,1)-size(E,1)+size(F,1); 

nRefScalar=(log(NF)-log(nF))/log(4);

if nRefScalar<0
    nRef=floor(nRefScalar);
    nRange=0:-1:nRef;
else
    nRef=ceil(nRefScalar);
    nRange=0:1:nRef;
end

nvR=nV.*ones(numel(nRange),1);
neR=nE.*ones(numel(nRange),1);

nfR=nF*4.^nRange';
for q=2:1:numel(nRange)
    if nRef>0                
        nvR(q)=nvR(q-1)+neR(q-1);        
    elseif nRef<0
        nvR(q)=(nvR(q-1)+X-nfR(q))/2;
    end    
    neR(q)=-X+nfR(q)+nvR(q);
end

% [nvR -neR nfR]
% sum([nvR -neR nfR],2)

nP=ceil(interp1(nfR(:),nvR(:),NF,'pchip'));