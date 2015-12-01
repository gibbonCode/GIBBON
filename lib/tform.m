function V2=tform(T,V)

%%
nDim_T=size(T,2);
nDim_V=size(V,2);

if nDim_T==(nDim_V+1); 
        VV=V;
        VV(:,end+1)=1;
        V2=(T*VV')';
        V2=V2(:,[1 2 3]);
elseif nDim_T==nDim_V;
        V2=(T*V')';
end

%%