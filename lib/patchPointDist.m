function [varargout]=patchPointDist(V,F,Vp,dL)

methodOpt='near-norm';

switch methodOpt
    case 'near-norm'
        
        [Nn,Vn]=patchNormal(F,V); %Get face normals and normal vector origins
        [Dn,indMin]=minDist(Vp,Vn); %Get distances to origins and find closest faces
        
        %Attempt to find orthogonal projection to closest face
        N=Nn(indMin,:); %Order the normals
        W=Vp-Vn(indMin,:);
        Dm=(dot(N,W,2));
        D=abs(Dm);        
        Dp=Dm(:,ones(1,3)).*N;
        Vpd=Vp-Dp;
        
        %Test if projection is valid
        TR = triangulation(F,V);
        B = cartesianToBarycentric(TR,indMin,Vpd);        
        logicOutside=any(B<-dL,2);
        
        %Replace invalid projections with face centre points
        D(logicOutside)=Dn(logicOutside); %Replace distances        
        Dp(logicOutside,:)=Vp(logicOutside,:)-Vn(indMin(logicOutside),:); %Replace difference vectors         
        
    otherwise
        error('Wrong method option chosen');
end

switch nargout
    case 1
        varargout{1}=D;
    case 2
        varargout{1}=D;
        varargout{2}=Dp;    
    case 3
        varargout{1}=D;
        varargout{2}=Dp;
        varargout{3}=logicOutside;
end