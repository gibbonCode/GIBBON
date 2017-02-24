function [F,V,C]=sweepLoft(V1,V2,n1,n2,Vg,numSteps,numTwist,plotOn)

%%

p1=Vg(1,:);
p2=Vg(end,:);


%%
if plotOn
    [Fq1,Vq1,Cq1]=quiver3Dpatch(p1(1),p1(2),p1(3),n1(1),n1(2),n1(3),[],[1 1]);
    [Fq2,Vq2,Cq2]=quiver3Dpatch(p2(1),p2(2),p2(3),n2(1),n2(2),n2(3),[],[1 1]);
    
    cFigure; hold on;
    xlabel('x');ylabel('y');zlabel('z');
    plotV(V1,'r-','lineWidth',3);
    plotV(V2,'b-','lineWidth',3);
    plotV(Vg,'g-','lineWidth',3);
    plotV(p1,'r.','MarkerSize',15);
    plotV(p2,'b.','MarkerSize',15);
    gpatch(Fq1,Vq1,'r');
    gpatch(Fq2,Vq2,'b');
    axis tight; axis equal; view(3);
    grid on; box on; grid on;
    drawnow;
end


%% Define allong curve coordinate systems

Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
U=Uf;
U(:,:,2)=Ub;
U=nanmean(U,3);
U=vecnormalize(U);

mean_V1=mean(V1,1);
V1m=vecnormalize(V1-mean_V1(ones(size(V1,1),1),:));
e3=n1;
e1=vecnormalize(V1m(1,:));
e2=vecnormalize(cross(e3,e1));
e1=vecnormalize(cross(e2,e3));
R1=[e1;e2;e3];

mean_V2=mean(V2,1);
V2m=vecnormalize(V2-mean_V2(ones(size(V2,1),1),:));
e3=n2;
e1=vecnormalize(V2m(1,:));
e2=vecnormalize(cross(e3,e1));
e1=vecnormalize(cross(e2,e3));
R2=[e1;e2;e3];

R_curve=repmat(eye(3,3),[1,1,size(Vg,1)]);
R_curve(:,:,1)=R1;
R=R1; 

if plotOn
    quiverTriad(Vg(1,:),R1,1); drawnow;
end

for q=2:size(Vg,1)
    a=U(q-1,:);
    b=U(q,:);
    
    theta=real(acos(dot(a,b))); %Complex if dot product is out of range [-1.1] due to precission issues
    w=vecnormalize(cross(b,a));
    [Rn]=vecAngle2Rot(theta,w);
    R=R*Rn;
    if plotOn
        quiverTriad(Vg(q,:),R,1);
        drawnow; 
    end
    R_curve(:,:,q)=R;
end

%% Create coordinate matrices

V1s=V1-p1(ones(size(V1,1),1),:);
V1s=V1s*R1';

V2s=V2-p2(ones(size(V2,1),1),:);
V2s=V2s*R2';

% Create coordinate matrices
X=linspacen(V1s(:,1),V2s(:,1),numSteps)';
Y=linspacen(V1s(:,2),V2s(:,2),numSteps)';
Z=linspacen(V1s(:,3),V2s(:,3),numSteps)';

Vs=[X(:) Y(:) Z(:)];

X=reshape(Vs(:,1),size(X));
Y=reshape(Vs(:,2),size(Y));
Z=reshape(Vs(:,3),size(Z));

%%

for q=1:1:numSteps

    V2p=[X(q,:)' Y(q,:)' Z(q,:)'];    
    V2p=(R_curve(:,:,q)'*V2p')';
    V2p=V2p+Vg(q*ones(size(V2p,1),1),:);    
    if plotOn==1
        plotV(V2p,'b-','LineWidth',1); drawnow; 
    end
    X(q,:)=V2p(:,1);
    Y(q,:)=V2p(:,2);
    Z(q,:)=V2p(:,3);
end

%%

V2p=[X(end,:)' Y(end,:)' Z(end,:)'];

[~,R]=rigidTransformationMatrixDirect(V2,V2p);

[theta,w]=rot2VecAngle(R);

%Flip angle
% if dot(U(end,:)',w)<0
%      theta=-theta;
% end

V2p_c=mean(V2p,1);
V2p=V2p-V2p_c(ones(size(V2p,1),1),:);
V2p=(R'*V2p')';
V2p=V2p+V2p_c(ones(size(V2p,1),1),:);

if plotOn==1
    plotV(V2p,'g--','lineWidth',3); drawnow;
end

%%

theta_w=linspace(0,1,size(Vg,1));
theta_step=theta*theta_w;

if numTwist>0
    theta_step_twist=linspace(0,2*pi*numTwist,size(Vg,1));
end

for q=1:1:size(Vg,1)    
    Vn=[X(q,:)' Y(q,:)' Z(q,:)'];
    Vn_mean=mean(Vn,1);
    Vn=Vn-Vn_mean(ones(size(Vn,1),1),:);
    
    [Rc]=vecAngle2Rot(theta_step(q),w);       
    Vn=(Rc'*Vn')';    
    
    if numTwist>0
        [Rc]=vecAngle2Rot(theta_step_twist(q),U(q,:));
        Vn=(Rc'*Vn')';
    end
    
    Vn=Vn+Vn_mean(ones(size(Vn,1),1),:);
    
    if plotOn
        plotV(Vn,'k-','lineWidth',1,'MarkerSize',25); drawnow; 
    end
    
    X(q,:)=Vn(:,1);
    Y(q,:)=Vn(:,2);
    Z(q,:)=Vn(:,3);
end


% %%
% 
% V1s=V1-p1(ones(size(V1,1),1),:);
% e3=n1;
% e1=vecnormalize(V1s(1,:));
% e2=vecnormalize(cross(e3,e1));
% e1=vecnormalize(cross(e2,e3));
% R1s=[e1;e2;e3];
% V1s=V1s*R1s';
% 
% V2m=mean(V2,1);
% V2s=V2-p2(ones(size(V2,1),1),:);
% e3=n2;
% e1=vecnormalize(V2s(1,:));
% e2=vecnormalize(cross(e3,e1));
% e1=vecnormalize(cross(e2,e3));
% R2s=[e1;e2;e3];
% V2s=V2s*R2s';
% 
% % Create coordinate matrices
% Xs=linspacen(V1s(:,1),V2s(:,1),numSteps)';
% Ys=linspacen(V1s(:,2),V2s(:,2),numSteps)';
% Zs=linspacen(V1s(:,3),V2s(:,3),numSteps)';
% 
% Vs=[Xs(:) Ys(:) Zs(:)];
% 
% Xs=reshape(Vs(:,1),size(Xs));
% Ys=reshape(Vs(:,2),size(Ys));
% Zs=reshape(Vs(:,3),size(Zs));
% 
% %% Define allong curve coordinate systems
% 
% Uf=[diff(Vg,1,1); nan(1,size(Vg,2))];
% Ub=-[nan(1,size(Vg,2)); flipud(diff(flipud(Vg),1,1))];
% U=Uf;
% U(:,:,2)=Ub;
% U=nanmean(U,3);
% U=vecnormalize(U);
% U(end,:)=n2; 
% 
% R=eye(3,3);
% e3=U(1,:);
% e2=[0 1 0];
% e1=vecnormalize(cross(e2,e3));
% e2=vecnormalize(cross(e3,e1));
% 
% R=[e1;e2;e3];
% 
% R_cell=cell(1,size(Vg,1));
% R1=R;
% R_cell{1}=R1s;
% for q=2:size(Vg,1)
%     a=U(q-1,:);
%     b=U(q,:);
%     
%     theta=real(acos(dot(a,b))); %Complex if dot product is out of range [-1.1] due to precission issues
%     w=vecnormalize(cross(b,a));
%     [Rn]=vecAngle2Rot(theta,w);
%     R=R*Rn;
%     if plotOn
%         quiverTriad(Vg(q,:),R,1);
%         drawnow; 
%     end
%     R_cell{q}=R;
% end
% 
% %% Fix for rotation
% 
% V1_now=[Xs(end,:)' Ys(end,:)' Zs(end,:)'];
% V2p=V1_now*R_cell{end};%[X(end,:)' Y(end,:)' Z(end,:)'];
% % V1_now=V1_now+p1(ones(size(V1_now,1),1),:);
% % V2p=V2p+p1(ones(size(V2p,1),1),:);
% V2p=V2p+Vg(size(Vg,1)*ones(size(V2p,1),1),:);
% 
% if plotOn
%     plotV(V2p,'r--','lineWidth',3);
%     drawnow; 
% end
% 
% [~,R]=rigidTransformationMatrixDirect(V2,V2p);
% 
% [theta,w]=rot2VecAngle(R);
% if dot(n2,w)<0
%     theta=-theta;
% end
% theta=theta+(2*pi)*numTwist;
% 
% % [R]=vecAngle2Rot(theta,w)
% 
% V2p_c=mean(V2p,1);
% V2p=((V2p-V2p_c(ones(size(V2p,1),1),:))*R)+V2p_c(ones(size(V2p,1),1),:);
% 
% % plotV(V2p,'g--','lineWidth',3);
% 
% %% Build coordinate matrices
% 
% X=zeros(size(Vg,1),size(V1,1));
% Y=X;
% Z=X;
% 
% theta_w=linspace(0,1,size(Vg,1));
% theta_step=theta*theta_w;
% Rc=eye(3,3);
% for q=1:1:size(Vg,1)
%     
%     V1_now=[Xs(q,:)' Ys(q,:)' Zs(q,:)'];    
%     Rn=R_cell{q};
%     w=Rn(3,:);
%     [Rc]=vecAngle2Rot(theta_step(q),w);        
%     Vpn=V1_now*Rn*Rc;
%     Vpn=Vpn+Vg(q*ones(size(Vpn,1),1),:);
%     if plotOn
%         plotV(Vpn,'k-','lineWidth',1,'MarkerSize',25);
%         drawnow; 
%     end
%     
%     X(q,:)=Vpn(:,1);
%     Y(q,:)=Vpn(:,2);
%     Z(q,:)=Vpn(:,3);
% end
% 
%% Override start and end with input curves

X(1,:)=V1(:,1);
Y(1,:)=V1(:,2);
Z(1,:)=V1(:,3);

X(end,:)=V2(:,1);
Y(end,:)=V2(:,2);
Z(end,:)=V2(:,3);

%% Create patch data

%Color data
c=(1:1:size(Z,1))';
C=c(:,ones(1,size(Z,2)));

%Create quad patch data
[F,V,C] = surf2patch(X,Y,Z,C);

%Close patch if required
I=[(2:size(Z,1))' (2:size(Z,1))' (1:size(Z,1)-1)' (1:size(Z,1)-1)'];
J=[ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) size(Z,2).*ones(size(Z,1)-1,1) ones(size(Z,1)-1,1)];
F_sub=sub2ind(size(Z),I,J);
F=[F;F_sub];

[C]=vertexToFaceMeasure(F,C);
C(end-size(F_sub,1):end,:)=C(end-size(F_sub,1):end,:)+0.5;
C=round(C);

if plotOn
    h=gpatch(F,V,C,'k',1);
    colormap gjet;
    camlight headlight;
    drawnow;
end

