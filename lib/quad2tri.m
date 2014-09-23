function [Ft,Vt]=quad2tri(Fq,Vq,triType)



switch triType
    case 'b' %backward slash
        Ft=[Fq(:,[1 2 3]);Fq(:,[3 4 1]);];
        Vt=Vq;
    case 'f' %Forward slash
        Ft=[Fq(:,[1 2 4]);Fq(:,[2 3 4]);];        
        Vt=Vq;
    case 'x' %Cross type
        Vm=zeros(size(Fq,1),size(Vq,2));
        for q=1:1:size(Vq,2)
            X=Vq(:,q);
            FX=X(Fq);
            if size(Fq,1)==1 %Treat special case of single face
                FX=FX';
            end
            Vm(:,q)=mean(FX,2);
        end
        %Join point sets
        Vt=[Vq;Vm];
        
        indVm=(size(Vq,1)+1):size(Vt,1);
        %Create faces
        Ft=[Fq(:,1) Fq(:,2) indVm(:);...
            Fq(:,2) Fq(:,3) indVm(:);...
            Fq(:,3) Fq(:,4) indVm(:);...
            Fq(:,4) Fq(:,1) indVm(:)];
end



