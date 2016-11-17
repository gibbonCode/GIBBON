function [Fc,Vc]=patchDetach(F,V,shrinkFactor)

Vc=zeros(size(F,1)*size(F,2),size(V,2));
for q=1:1:size(V,2)
    X=V(:,q);
    if size(F,1)==1
        FX=X(F)';
    else
        FX=X(F);
    end
    FX_mean=mean(FX,2);
    FX=((FX-FX_mean)*shrinkFactor)+FX_mean;
    Vc(:,q)=FX(:);
end
    
Fc=reshape(1:size(Vc,1),size(F,1),size(F,2));