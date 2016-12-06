function [Q]=pointCloudPrinDir(V)

[U_svd,~,~]=svd(V',0);
%U_svd=U_svd./norm(U_svd);
Q=U_svd';