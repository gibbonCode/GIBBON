function [Fn,Vn]=subTriCentre(F,V,L)

% function function [Fn,Vn]=subTriCentre(F,V,L)
% ------------------------------------------------------------------------
% This function splits the faces defined by L up into three by introducing
% a central node. 
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/06/03
%------------------------------------------------------------------------

numPoints=size(V,1);
Fs=F(L,:); %Faces to split

%Get face centres
X=V(:,1); Y=V(:,2); Z=V(:,3);
if nnz(L)==1 %treat different behaviour for since face problem
    XF=mean(X(Fs),1);
    YF=mean(Y(Fs),1);
    ZF=mean(Z(Fs),1);
else
    XF=mean(X(Fs),2);
    YF=mean(Y(Fs),2);
    ZF=mean(Z(Fs),2);
end
V_add=[XF(:) YF(:) ZF(:)];
numPointsAdd=size(V_add,1);

Vn=[V;V_add];
indNew=(numPoints+1:numPoints+numPointsAdd)';
F_add=[Fs(:,2) indNew Fs(:,1);...
       Fs(:,3) indNew Fs(:,2);...
       Fs(:,1) indNew Fs(:,3)];

Fn=[F(~L,:);F_add];