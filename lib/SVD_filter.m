function [Zm]=SVD_filter(Z,P,T)

%Computer SVD
[U,S,V] = svd(Z);

%Normalising singular values
Ss=diag(S)-min(diag(S));
Ss=Ss./max(Ss);

%Creating smoothening parameters
p_max=P(1); %1=No blurring
p_min=P(2); %0=straight line fit
p=(Ss.*(p_max-p_min))+p_min; %Scale towards singular values

Zm=nan(size(Z));
for i=1:1:size(U,2);
    v=V(:,i); u=U(:,i); s=S(i,i); %components    
    if Ss(i)<T %Filter after threshold
                us = csaps(1:numel(u),u,p(i),1:numel(u))'; %Smooth u
                vs = csaps(1:numel(v),v,p(i),1:numel(v))'; %Smooth v
    else
        vs=v; us=u; %Keep unsmoothened
    end
    z=us*s*vs'; %sub-data
    Zm(:,:,i)=z;
end
Zm=sum(Zm,3);
