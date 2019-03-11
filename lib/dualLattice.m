function [FT,VT,CT]=dualLattice(E,V,shrinkFactor,cladOpt)

% WORK IN PROGRESS

%%

[Ec,Vc]=patchDetach(E,V,shrinkFactor);

[F]=element2patch(E); % [faceMat]=patchEdges(F,0);
[F_Ec]=element2patch(Ec); % Ec=patchEdges(Fc,0);

[Fs,indSort]=sort(F,2);
[~,indUni1,indUni2,uniCount,IND_MAP]=unique_map(Fs,'rows');
IND_MAP=sort(IND_MAP,1,'descend');
IND_MAP=full(IND_MAP(1:2,:))';

I=(1:1:size(F,1))';
INDSort=sub2ind(size(F),I(:,ones(size(indSort,2),1)),indSort);

indSort2=indSort;
indSort2=indSort2(indUni1,:); %Unique
indSort2=indSort2(indUni2,:); %Expand

[~,indUnsort]=sort(indSort2,2);
I=(1:1:size(F,1))';
INDUnsort=sub2ind(size(F),I(:,ones(size(indUnsort,2),1)),indUnsort);

F_Ec=F_Ec(INDSort);
F_Ec=F_Ec(INDUnsort);

F=F(indUni1,:);
if numel(shrinkFactor)>1
    if numel(shrinkFactor)==size(E,1) %Shrink factor specified on faces
        [shrinkFactor_V]=faceToVertexMeasure(E,V,shrinkFactor); %Convert to nodal metric
    else
        shrinkFactor_V=shrinkFactor;
    end
    shrinkFactor_E=mean(shrinkFactor_V(F),2); %Convert to edge metric
else
    shrinkFactor_E=shrinkFactor;
end

Fq=NaN(size(IND_MAP,1),6);
Fq(uniCount==2,:)=[F_Ec(IND_MAP(uniCount==2,1),:) F_Ec(IND_MAP(uniCount==2,2),:)];

if any(uniCount==1)
    [Fcc,Vcc]=patchDetach(F,V,shrinkFactor_E);
    Fcc=Fcc(indUni2,:); %Expand
    Vq=[Vc;Vcc];
    F1=fliplr(F_Ec(IND_MAP(uniCount==1,1),:));
    F2=fliplr(Fcc(IND_MAP(uniCount==1,1),:))+size(Vc,1);
    Fq(uniCount==1,:)=[F1 F2];
else
    Vq=Vc;
end

Fq=[Fq(:,[1 2]) Fq(:,[5 4]);...
    Fq(:,[2 3]) Fq(:,[6 5]);...
    Fq(:,[3 1]) Fq(:,[4 6])];

Fq=fliplr(Fq);

if cladOpt==1
    F=F(indUni2,:); %Expand
    [Fq2,Vq2,Fc,~]=dualClad(F(IND_MAP(uniCount==1,1),:),V,shrinkFactor,1);
    
    % [~,~,Nv]=patchNormal(Fc,Vq2);
    
    D=patchEdgeLengths(Fc,Vq2);
    D=reshape(D',size(Fc,2),size(Fc,1))';
    D=min(D,[],2);
    D=faceToVertexMeasure(Fc,Vq2,D);
    
    [~,indUni_F2,~]=unique(F2(:));
    distEdge=sqrt(sum((Vq(F1(indUni_F2),:)-Vq(F2(indUni_F2),:)).^2,2));
    q=((distEdge-D)./distEdge);
    Vq(F2(indUni_F2),:)=((1-q).*Vq(F1(indUni_F2),:))+((q).*Vq(F2(indUni_F2),:));
    
    Vq2i=Vq2;
    F2=fliplr(F2);
    Vq2i(Fc(indUni_F2),:)=Vq(F2(indUni_F2),:);
    
    Vq3=[Vq2;Vq2i];
    
    Fq2_p=fliplr(Fq2)+size(Vq2,1);
    
    Fq2_2=[Fq2(:,[3 2]) Fq2(:,[2 3])+size(Vq2,1)];%
    Fq2_3=[Fq2(:,[1 4]) Fq2(:,[4 1])+size(Vq2,1)];%
    
    Fq2_sides=[Fq2_p; Fq2_2; Fq2_3];
    
    Fq3=[Fq2_sides; Fq2;];
    Cq3=[2*ones(size(Fq2_sides,1),1); 3*ones(size(Fq2,1),1);];
    
    Cc=4*ones(size(Fc,1),1);
    
    FT=[Fq(:,[1 2 3]); Fq(:,[3 4 1]); Fq3(:,[1 2 3])+size(Vq,1); Fq3(:,[3 4 1])+size(Vq,1); Fc+size(Vq,1)];
    VT=[Vq; Vq3];
    CT=[ones(size(Fq,1)*2,1); repmat(Cq3,2,1); Cc; ];        
    [FT,VT]=mergeVertices(FT,VT,5);
else
   FT=[Fq(:,[1 2 3]); Fq(:,[3 4 1]); fliplr(F2)];   
   VT=Vq;
   CT=[ones(size(Fq,1)*2,1);2*ones(size(F2,1),1)];
end



%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
