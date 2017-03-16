function [varargout]=polyNormal(V_poly)

% [N,Vn,Nv]=polyNormal(V_poly)
% ------------------------------------------------------------------------
% Normals are derived based on cross product of triangle edge vectors. Each
% triangle is constructed using the two points of its face and the mean of
% the face. 
%
%
% To do: Check for co-linear edges ?
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2014/06/02 %Updated general lay-out and commenting
% 2015/09/22 %Fixed to allow for 2D patch data
% 2016/11/15 %Added handling of non-triangular faces
%------------------------------------------------------------------------

%%

zDir=[0 0 1];

if size(V_poly,2)==2 %Cope with 2D
    V_poly(:,3)=0;
end

[Vn,U]=polyDerivative(V_poly,1);
N=cross(zDir(ones(size(U,1),1),:),U);
N=vecnormalize(N);

varargout{1}=N(:,[1 2]); 
varargout{2}=Vn(:,[1 2]); 

if nargout==3
    [~,U]=polyDerivative(V_poly,3);
    Nv=cross(zDir(ones(size(U,1),1),:),U);
    Nv=vecnormalize(Nv);
    varargout{3}=Nv(:,[1 2]);
end

end

function [V,U]=polyDerivative(Vg,dirOpt)

switch dirOpt
    case 1 %Forward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=diff(Vg,1,1);        
    case 2 %Backward
        V=(Vg(1:end-1,:)+Vg(2:end,:))/2;
        U=flipud(diff(flipud(Vg),1,1));        
    case 3 %Central
        V=Vg;
        Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
        Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
        U=Uf;
        U(:,:,2)=Ub;
        U=nanmean(U,3);        
end
U=vecnormalize(U);

end