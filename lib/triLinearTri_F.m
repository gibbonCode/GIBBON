function F_cell=triLinearTri_F(VX,Vx,TRI)

F_cell=cell(size(TRI,1),1);
for qc=1:1:numel(F_cell);
    indNow=TRI(qc,:);
    X=VX(indNow,:); %Initial coordinates
    x=Vx(indNow,:); %Current coordinates
    F=triLinearTri_subF(X,x); %The deformation gradient tensor
    F_cell{qc}=F; %Store in cell array
end

    
end

function F=triLinearTri_subF(X,x)
%% Define shape function set

% syms r s;
% N=[1-r-s; r; s;];

%% COMPUTE DERIVATIVES

% %Compute derivatives of shape functions with respect to RST
% dN_dRST=zeros(size(N,1),size(X,2));
% for q=1:1:size(N,1);
%     dN_dRST(q,:)=[diff(N(q),r) diff(N(q),s) diff(N(q),t)];
% end

%Hardcoded derivative
dN_dRST =[-1    -1 ;...
           1     0 ;...
           0     1 ;];

%Compute derivatives of initial position vectors with respect to shape functions
dX_dRST=zeros(3,3);
for q=1:1:4;
    dX_dRST=dX_dRST+(X(q,:)'*dN_dRST(q,:));
end
dX_dRST_invT=inv(dX_dRST).';

%Compute derivatives of current position vectors with respect to shape functions
dx_dRST=zeros(3,3);
for q=1:1:4;
    dx_dRST=dx_dRST+(x(q,:)'*dN_dRST(q,:));
end
dx_dRST_invT=inv(dx_dRST)';

%Compute derivatives of shape functions with respect to initial position vectors
dN_dX=zeros(4,3);
for q=1:1:4;
    dN_dX(q,:)=(dX_dRST_invT*dN_dRST(q,:)')';
%     dN_dX(q,:)=(dX_dRST'\dN_dRST(q,:)')';
end

%Compute derivatives of shape functions with respect to current position vectors
dN_dx=zeros(4,3);
for q=1:1:4;
    dN_dx(q,:)=(dx_dRST_invT*dN_dRST(q,:)')';
end

%% DERIVE THE DEFORMATION GRADIENT TENSOR
F=zeros(3,3);
for q=1:1:size(x,1)
    F=F+(x(q,:)'*dN_dX(q,:));
end

end
