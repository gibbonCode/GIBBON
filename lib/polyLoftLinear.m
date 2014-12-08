function [F,V]=polyLoftLinear(Vc_start,Vc_end,cPar)

X=linspacen(Vc_start(:,1),Vc_end(:,1),cPar.numSteps)';
Y=linspacen(Vc_start(:,2),Vc_end(:,2),cPar.numSteps)';
Z=linspacen(Vc_start(:,3),Vc_end(:,3),cPar.numSteps)';
c=(1:1:size(Z,1))';
C=c(:,ones(1,size(Z,2)));

%Create quad patch data
[F,V,C] = surf2patch(X,Y,Z,C);

%Close patch if required
if cPar.closeLoopOpt
    I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
    J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
    F_sub=sub2ind(size(Z),I,J);
    F=[F;F_sub];
    [C]=vertexToFaceMeasure(F,C);
    C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5; 
else
    [C]=vertexToFaceMeasure(F,C);
end
C=round(C);

switch cPar.patchType
    case 'quad' 
    case 'tri_slash' %Convert quads to triangles by slashing
        F=[F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)];
    case 'tri' %Convert quads to approximate equilateral triangles

        logicSlashType=repmat(iseven(C),2,1);
        
        Xi=X;
        x=X(1,:);
        dx=diff(x)/2;
        dx(end+1)=(x(1)-x(end))/2;
        for q=2:2:size(X,1)
            X(q,:)=X(q,:)+dx;
        end
        if ~cPar.closeLoopOpt
            X(:,1)=Xi(:,1);
            X(:,end)=Xi(:,end);
        end
        
        Yi=Y;
        y=Y(1,:);
        dy=diff(y)/2;
        dy(end+1)=(y(1)-y(end))/2;
        for q=2:2:size(Y,1)
            Y(q,:)=Y(q,:)+dy;
        end
        if ~cPar.closeLoopOpt
            Y(:,1)=Yi(:,1);
            Y(:,end)=Yi(:,end);
        end
        
        Zi=Z;
        z=Z(1,:);
        dz=diff(z)/2;
        dz(end+1)=(z(1)-z(end))/2;
        for q=2:2:size(Z,1)
            Z(q,:)=Z(q,:)+dz;
        end
        if ~cPar.closeLoopOpt
            Z(:,1)=Zi(:,1);
            Z(:,end)=Zi(:,end);
        end
        
        V=[X(:) Y(:) Z(:)];
        
        F1=[F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)]; 
        F2=[F(:,1) F(:,4) F(:,2); F(:,2) F(:,4) F(:,3)];
        F=[F1(~logicSlashType,:);F2(logicSlashType,:)];
        
%         C=repmat(C,2,1);
%         C=[C(~logicSlashType,:);C(logicSlashType,:)];
        
end

