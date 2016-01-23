function [F,V]=patchcylinder(r,nr,h,nz,ptype)

% [THETA,Z] = meshgrid(linspace(0,2*pi,nr+1),h.*linspace(-0.5,0.5,nz));
% THETA=THETA(:,1:end-1); Z=Z(:,1:end-1);
% [X,Y,Z] = pol2cart(THETA,r.*ones(size(Z)),Z);
% [F,V] = surf2patch(X,Y,Z);
% I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
% J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
% F_sub=sub2ind(size(Z),I,J);
% F=[F;F_sub];
% 
% switch ptype
%     case 'quad' %Meshed cylinder
% 
%     case 'tri' %Triangulated cylinder
%         F=fliplr([F(:,1) F(:,3) F(:,2); F(:,1) F(:,4) F(:,3)]);
%     otherwise
%         warning('wrong input for argument ptype, valid inputs are quad and tri');
% end

t=linspace(0,2*pi,nr+1);
t=t(1:end-1);
x=r*cos(t);
y=r*sin(t);
Vc=[x(:) y(:)];
Vc(:,3)=0; 

cPar.numSteps=nz;
cPar.depth=h; 
cPar.patchType=ptype; 
cPar.dir=0;
cPar.closeLoopOpt=1; 

[F,V]=polyExtrude(Vc,cPar);

end