function [Fq,Vq,Fc,Vc]=dualClad(F,V,shrinkFactor,cladMethod)

[Fc,Vc]=patchDetach(F,V,shrinkFactor);

[E]=patchEdges(F,0);
Es=sort(E,2);
[~,indUni1,indUni2,uniCount,IND_MAP]=unique_map(Es,'rows');
IND_MAP=sort(IND_MAP,1,'descend');
IND_MAP=full(IND_MAP(1:2,:))';

switch cladMethod
    case 1
        E=E(indUni1,:);
        if numel(shrinkFactor)>1
            if numel(shrinkFactor)==size(F,1) %Shrink factor specified on faces
                [shrinkFactor_V]=faceToVertexMeasure(F,V,shrinkFactor); %Convert to nodal metric
            else
                shrinkFactor_V=shrinkFactor;
            end
            shrinkFactor_E=mean(shrinkFactor_V(E),2); %Convert to edge metric
        else
            shrinkFactor_E=shrinkFactor;
        end
        Ec=patchEdges(Fc,0);
        Fq=zeros(size(IND_MAP,1),4);
        Fq(uniCount==2,:)=fliplr([Ec(IND_MAP(uniCount==2,1),:) Ec(IND_MAP(uniCount==2,2),:)]);
        if any(uniCount==1)
            [Ecc,Vcc]=patchDetach(E,V,shrinkFactor_E);
            Ecc=Ecc(indUni2,:); %Expand
            Ec=patchEdges(Fc,0);
            Vq=[Vc;Vcc];
            Fq(uniCount==1,:)=[Ec(IND_MAP(uniCount==1,1),[2 1]) Ecc(IND_MAP(uniCount==1,1),:)+size(Vc,1);];
        else
            Vq=Vc;
        end
    case 2
        [E_uni]=patchEdges(F,1);
        if numel(shrinkFactor)>1
            if numel(shrinkFactor)==size(F,1) %Shrink factor specified on faces
                [shrinkFactor_V]=faceToVertexMeasure(F,V,shrinkFactor); %Convert to nodal metric
            else
                shrinkFactor_V=shrinkFactor;
            end
            shrinkFactor_E=mean(shrinkFactor_V(E_uni),2); %Convert to edge metric
        else
            shrinkFactor_E=shrinkFactor;
        end
        logicFlip=~all(E_uni(indUni2,:)==E,2);
        [Ecc,Vcc]=patchDetach(E_uni,V,shrinkFactor_E);
        Ecc=Ecc(indUni2,:); %Expand
        Ecc(logicFlip,:)=fliplr(Ecc(logicFlip,:));
        Ec=patchEdges(Fc,0);
        Fq=fliplr([Ecc(:,[2 1])+size(Vc,1) Ec ]);
        Vq=[Vc;Vcc];
    case 3
        
        %%%%%%% Like method 2
        [E_uni]=patchEdges(F,1);
        if numel(shrinkFactor)>1
            if numel(shrinkFactor)==size(F,1) %Shrink factor specified on faces
                [shrinkFactor_V]=faceToVertexMeasure(F,V,shrinkFactor); %Convert to nodal metric
            else
                shrinkFactor_V=shrinkFactor;
            end
            shrinkFactor_E=mean(shrinkFactor_V(E_uni),2); %Convert to edge metric
        else
            shrinkFactor_E=shrinkFactor;
        end
        logicFlip=~all(E_uni(indUni2,:)==E,2);
        [Ecc,Vcc]=patchDetach(E_uni,V,shrinkFactor_E);
        Ecc=Ecc(indUni2,:); %Expand
        Ecc(logicFlip,:)=fliplr(Ecc(logicFlip,:));
        Ec=patchEdges(Fc,0);
        Fq=fliplr([Ecc(:,[2 1])+size(Vc,1) Ec ]);
        Vq=[Vc;Vcc];
        %%%%%
        
        Fqc=[Ec(IND_MAP(uniCount==2,1),:) Ec(IND_MAP(uniCount==2,2),:)];
        
        e1=[Fqc(:,1) Fqc(:,4)];
        e2=[Fqc(:,2) Fqc(:,3)];
        
        a=Vc(e1(:,2),:)-Vc(e1(:,1),:);
        A=Vc(e1(:,1),:);
        b=V(E(IND_MAP(uniCount==2,1),2),:)-V(E(IND_MAP(uniCount==2,1),1),:);
        B=V(E(IND_MAP(uniCount==2,1),1),:);
        [~,X1]=vecPairClosestPoint(a,b,A,B,0);
        
        a=Vc(e2(:,2),:)-Vc(e2(:,1),:);
        A=Vc(e2(:,1),:);
        
        [~,X2]=vecPairClosestPoint(a,b,A,B,0);
        
        ind=Ecc(IND_MAP(uniCount==2,1),:);
        VX=[X1;X2];
        Vq(ind+size(Vc,1),:)=VX;
end
