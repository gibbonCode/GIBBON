function [varargout]=patch_dual(V,F)

% function [Vd,Fd,Fds]=patch_dual(V,F)
% ------------------------------------------------------------------------
% Computes the dual of the input tesselation defined by the vertices V and
% faces F. 
%
%
% 
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
% 2013/04/28
% 2015/07/07 Updated for GIBBON
%
% TO DO allow multiple iterations. Add "book keeping" of color/index data
%------------------------------------------------------------------------
%%

%Get input tesselation vertex-normals
[~,~,Nv]=patchNormal(F,V);

%Create patch indices
[IND_F]=tesIND(F,V,0);

%Face centre point coordinates
X=V(:,1); Y=V(:,2); Z=V(:,3);
Vd=[mean(X(F),2) mean(Y(F),2) mean(Z(F),2)];

%Creating arrays for faces
[I,J,v] = find(IND_F);
Xfd=accumarray({I,J},Vd(v,1),size(IND_F),[],NaN); 
Yfd=accumarray({I,J},Vd(v,2),size(IND_F),[],NaN); 
Zfd=accumarray({I,J},Vd(v,3),size(IND_F),[],NaN); 

Xfd_mean=nanmean(Xfd,2); 
Yfd_mean=nanmean(Yfd,2); 
Zfd_mean=nanmean(Zfd,2); 

Xfd=Xfd-Xfd_mean(:,ones(1,(size(Xfd,2))));
Yfd=Yfd-Yfd_mean(:,ones(1,(size(Yfd,2))));
Zfd=Zfd-Zfd_mean(:,ones(1,(size(Zfd,2))));

%Determine face order
Xfdn=Xfd; Yfdn=Yfd; Zfdn=Zfd;
for q=1:1:size(Xfd,1)
    Vc=[Xfd(q,:); Yfd(q,:); Zfd(q,:)];
    v1=[Xfd(q,1),Yfd(q,1),Zfd(q,1)]-[Xfd(q,2),Yfd(q,2),Zfd(q,2)]; [v1]=vecnormalize(v1);
    v2=[Xfd(q,2),Yfd(q,2),Zfd(q,2)]-[Xfd(q,3),Yfd(q,3),Zfd(q,3)]; [v2]=vecnormalize(v2);
    v3=cross(v1,v2); [v3]=vecnormalize(v3);
    v2=cross(v1,v3); [v2]=vecnormalize(v2);
    DCM=[v1; v2; v3];
    Vcn=(DCM*Vc)';
    Xfdn(q,:)=Vcn(:,1);
    Yfdn(q,:)=Vcn(:,2);
    Zfdn(q,:)=Vcn(:,3);
end

%Fix face order
theta_n=atan2(Yfdn,Xfdn); %[theta_n,~,~] = cart2pol(Xfdn,Yfdn,Zfdn);
[~,J_sort]=sort(theta_n,2);
% I_sort=(1:1:size(J_sort,1))'*ones(1,size(J_sort,2));
I_sort=(1:1:size(J_sort,1))'; 
I_sort=I_sort(:,ones(1,size(J_sort,2)));
IND_sort = sub2ind(size(J_sort),I_sort,J_sort);

% Creating faces matrix
Fds=IND_F(IND_sort);

%Splitting up into seperate face types and fix face normals
n_sum=sum(Fds>0,2);
face_num_types=unique(n_sum);
face_num_types=face_num_types(face_num_types>0);
Fd=cell(1,numel(face_num_types));
for q=1:1:numel(face_num_types)
    %Get faces
    Lf=(n_sum==face_num_types(q)); %logic for current faces    
    F_now=Fds(Lf,1:face_num_types(q)); %The current face set
    
    %Flip face orientation if required
    Nv_now=Nv(Lf,:); %Appropriate face normals based on input mesh    
    [N_now]=patchNormal(F_now,Vd); %Current face normals    
    D=sqrt(sum((N_now-Nv_now).^2,2)); %Difference metric    
    logicFlip=D>1;%max(eps(N_now(:))); %Logic for faces that require flipping    
    F_now(logicFlip,:)=fliplr(F_now(logicFlip,:)); %Flip faces
    
    ind_f=find(Lf); 
    indFlip=ind_f(logicFlip);
    Fds(indFlip,1:face_num_types(q))=fliplr(Fds(indFlip,1:face_num_types(q))); 
    
    %Store faces in cell array    
    Fd{q}=F_now; 
end

varargout{1}=Vd;
varargout{2}=Fd;
varargout{3}=Fds;
