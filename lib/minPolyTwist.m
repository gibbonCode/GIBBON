function [varargout]=minPolyTwist(V1,V2)

%Get number of points for curves
n1=size(V1,1);
n2=size(V2,1);

if n1==n2
    
    %Fix overall direction if required
    [L1]=isPolyClockwise(V1);
    [L2]=isPolyClockwise(V2);
    flipReq=L1~=L2;
    
    if flipReq
        flipFlag=1;
        V2t=flipud(V2);
    else
        flipFlag=0;
        V2t=V2;
    end
    
    %Find point order for minimum torsion
    SSQD=zeros(n1,1);
    for q=1:1:n1
        if q>1
            indSort=[q:n1 1:q-1];
        else
            indSort=1:n1;
        end
        V2s=V2t(indSort,:);
        SSQD(q)=sum(sqrt(sum((V1-V2s).^2,2)).^2);
    end
    
    %Get sort order
    [~,qMin]=min(SSQD);
    if qMin>1
        indSort=[qMin:n1 1:qMin-1];
    else
        indSort=1:n1;
    end
    
    %Fix order if flipping was required
    if flipReq
        indSort=n1-(indSort-1);
    end
    
    %Reorder points
    V2f=V2(indSort,:);
    
    switch nargout
        case 1
            varargout{1}=V2f;
        case 2
            varargout{1}=V2f;
            varargout{2}=indSort;
        case 3
            varargout{1}=V2f;
            varargout{2}=indSort;
            varargout{3}=flipFlag;            
    end
else
    error('size(V1,1) should equal size(V2,1)');
end

