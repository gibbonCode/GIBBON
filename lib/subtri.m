function [Fs,Vs]=subtri(F,V,n,uniqueOpt)

% function [Fs,Vs]=subtri(F,V,n,uniqueOpt)
% ------------------------------------------------------------------------
% Sub-triangulates the triangles defined by the patch format data F (faces)
% and V (vertices). Creates n addition points on the edges of the initial
% triangles, thus it creates (n+1).^2 triangles per original triangle.
% Two methods are implemented one which iteratively splits edges and can
% thus only be used if log2(n+1) is a positive integer. This method
% guarantees that all points are unique. If log2(n+1) is not a positive
% integer points are seeded on all edges initially creating non-unique
% points. However if uniqueOpt is 1 (default if not provided) these points
% are suppressed using the unique command. 
% 
%
% %% EXAMPLE 
% [X,Y] = meshgrid(linspace(-10,10,15)); 
% Z = sinc(sqrt((X/pi).^2+(Y/pi).^2)); 
% F = delaunay(X,Y); V=[X(:) Y(:) Z(:)]; C=mean(Z(F),2);
% 
% n=2; 
% [Fs,Vs]=subtri(F,V,n); 
% Vs(:,3)=sinc(sqrt((Vs(:,1)/pi).^2+(Vs(:,2)/pi).^2)); Z=Vs(:,3);Cs=mean(Z(Fs),2);
% 
% figure('units','normalized','Position',[0 0 1 1],'Color','w'); colordef('white'); 
% subplot(1,2,1);patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'FaceAlpha',0.5,'EdgeColor','k','LineWidth',2); hold on; 
% axis tight; axis square; grid on; hold on; view(3); axis off; 
% title('Original','FontSize',20); 
% subplot(1,2,2);patch('Faces',Fs,'Vertices',Vs,'FaceColor','flat','CData',Cs,'FaceAlpha',0.5,'EdgeColor','k','LineWidth',0.5); hold on; 
% axis tight; axis square; grid on; hold on; view(3); axis off; 
% title(['n=',num2str(n)],'FontSize',20);
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
%
% 2010/06/01 Created
% 2014/02/27 Added splitting method
% ------------------------------------------------------------------------

if nargin==3
    uniqueOpt=1;
end

%%

%Check if n can be achieved through splitting
nSplitIterations=log2(n+1); %Check for integer solution
logicInteger=abs(round(nSplitIterations)-nSplitIterations)<eps(nSplitIterations);

if uniqueOpt==1 && logicInteger
    subMethod='split';
else
    subMethod='seed';
end

switch subMethod
    case 'split' %iteratively split edges (no double points created, unique operation avoided)
            Fs=F; Vs=V;
        for q=1:1:nSplitIterations
            F=Fs; V=Vs;
            E=[F(:,[1 2]); F(:,[2 3]);  F(:,[3 1])]; %Edges matrix
            Es=sort(E,2); %Sorted edges matrix
            [~,ind1,~]=unique(Es,'rows');
            E=E(ind1,:);
            
            numPoints = size(V,1);
            numEdges = size(E,1);
            
            % Get indices of the three edges associated with each face            
            A = sparse(E(:,1),E(:,2),(1:numEdges)+numPoints,numPoints,numPoints,numEdges);            
            A = max(A,A'); %Copy symmetric
            
            %Indices for A matrix
            indA_12=F(:,1)+(F(:,2)-1)*numPoints;
            indA_23=F(:,2)+(F(:,3)-1)*numPoints;
            indA_31=F(:,3)+(F(:,1)-1)*numPoints;
            
            %Get indices for vertex array
            indV_12=full(A(indA_12));
            indV_23=full(A(indA_23));
            indV_31=full(A(indA_31));
            
            %Create faces array
            Fs=[[F(:,1)  indV_12 indV_31];...
                [F(:,2)  indV_23 indV_12];...
                [F(:,3)  indV_31 indV_23];...
                [indV_12 indV_23 indV_31]];
            
            %Create vertex array
            Vn=0.5*(V(E(:,1),:)+V(E(:,2),:)); %new mid-edge points
            Vs = [V; Vn]; %Join point sets
        end
    case 'seed' %Seed edge points and remove doubles (more memory intensive)
        nDim=size(V,2);
        if nDim==2
            V(:,3)=0;
        end
        
        no_faces=size(F,1);
        
        EF=edges(triangulation(F,V));
        X=V(:,1); Y=V(:,2); Z=V(:,3);
        X_EF=X(EF); Y_EF=Y(EF); Z_EF=Z(EF);
        edgeLengthMin=min(sqrt(sum(([X_EF(:,1) Y_EF(:,1) Z_EF(:,1)]-[X_EF(:,2) Y_EF(:,2) Z_EF(:,2)]).^2,2)));
        if edgeLengthMin>0
            scaleFactor=1./edgeLengthMin;
        else
            scaleFactor=1;
        end
        V=V.*scaleFactor;
        
        X=V(:,1); Y=V(:,2); Z=V(:,3);
        
        %Creating "Edge indices"
        vsize=size(V,1).*ones(1,3); %Arbitrary
        IJ_1=[F(:,1) F(:,2)];
        edge_IND1= sub2ind(vsize,IJ_1(:,1),IJ_1(:,2));
        IJ_2=[F(:,2) F(:,3)];
        edge_IND2= sub2ind(vsize,IJ_2(:,1),IJ_2(:,2));
        IJ_3=[F(:,3) F(:,1)];
        edge_IND3= sub2ind(vsize,IJ_3(:,1),IJ_3(:,2));
        
        E=[edge_IND1 edge_IND2 edge_IND3]; %Edges matrix
        [E_unique,~,ind_uni_2]=unique(E); %To remove doubles
        IND_E_uni=reshape(ind_uni_2,size(E));
        [I,J] = ind2sub(vsize,E_unique);
        
        % Creating points on edges
        A=X(I); B=X(J); A=A(:); B=B(:);
        Xe=(A*ones(1,n+2))+((B-A)./(n+1))*(0:1:n+1);
        Xe=Xe(:,1:end);
        
        A=Y(I); B=Y(J); A=A(:); B=B(:);
        Ye=(A*ones(1,n+2))+((B-A)./(n+1))*(0:1:n+1);
        Ye=Ye(:,1:end);
        
        A=Z(I); B=Z(J); A=A(:); B=B(:);
        Ze=(A*ones(1,n+2))+((B-A)./(n+1))*(0:1:n+1);
        Ze=Ze(:,1:end);
        
        %Creating points on faces
        Xi=[]; Yi=[]; Zi=[];
        for i=0:1:n+1
            if i==0
                X1=Xe(IND_E_uni(:,1),end); Y1=Ye(IND_E_uni(:,1),end); Z1=Ze(IND_E_uni(:,1),end);
            else
                no_steps=i+1;
                pfrom=Xe(IND_E_uni(:,1),end-i); pto=Xe(IND_E_uni(:,2),i+1);
                X1=(pfrom*ones(1,no_steps))+((pto-pfrom)./(no_steps-1))*(0:1:no_steps-1);
                pfrom=Ye(IND_E_uni(:,1),end-i); pto=Ye(IND_E_uni(:,2),i+1);
                Y1=(pfrom*ones(1,no_steps))+((pto-pfrom)./(no_steps-1))*(0:1:no_steps-1);
                pfrom=Ze(IND_E_uni(:,1),end-i); pto=Ze(IND_E_uni(:,2),i+1);
                Z1=(pfrom*ones(1,no_steps))+((pto-pfrom)./(no_steps-1))*(0:1:no_steps-1);
            end
            X1=X1(:,1:end);Y1=Y1(:,1:end);  Z1=Z1(:,1:end);
            Xi=[Xi; X1(:)]; Yi=[Yi; Y1(:)]; Zi=[Zi; Z1(:)];
        end
        
        %Setting up new faces matrix
        no_points_per_face=(0.5.*n+1).*(n+3);
        IND_F1=1+(no_faces.*((1:1:no_points_per_face)-1));
        IND_IND_F1=1:1:length(IND_F1);
        line_pairs=[IND_IND_F1(2:end-1)' IND_IND_F1(3:end)'];
        N=1:n+2;
        IND_skip_edge=(N.*(N+1))./2;
        L_keep=~ismember(line_pairs(:,1),IND_skip_edge);
        line_pairs=line_pairs(L_keep,:);
        F11=[line_pairs(:,1) (1:1:size(line_pairs,1))' line_pairs(:,2)];
        line_pairs=line_pairs(1:end-(n+1),:);
        F12=[line_pairs sum(line_pairs,2)-(0:1:size(line_pairs,1)-1)'];
        F1=[IND_F1(F11); IND_F1(F12)];
        A=ones(size(F1,1),1)*(0:1:size(F,1)-1);A=A(:)*ones(1,3);
        Fs=repmat(F1,[size(F,1),1])+A;
        
        Vs=[Xi(:) Yi(:) Zi(:)];
        
        if uniqueOpt 
            %this is based on rounding to 6th signigicant digit to avoid this use the split method            
            numKeep=5; %Number of significant digits to keep
            [~,ind_uni_1,ind_uni_2]=unique(pround(Vs,numKeep),'rows');
%             try
%                 [~,ind_uni_1,ind_uni_2]=unique(round(Vs,numKeep,'significant'),'rows');
%             catch
%                 [~,ind_uni_1,ind_uni_2]=unique(sround(Vs,numKeep),'rows');
%             end
            Vs=Vs(ind_uni_1,:);
            Fs=ind_uni_2(Fs);
        end
        Vs=Vs./scaleFactor;
        if nDim==2
           Vs=Vs(:,[1 2]); 
        end
end

end