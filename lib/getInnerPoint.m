function [V_in]=getInnerPoint(F,V,searchRadius,voxelSize,plotOn)

%% Get an inner point
%Convert patch data to image
[M,G,~]=triSurf2Im(F,V,voxelSize);
imOrigin=G.origin; %Image origin can be used to allign image with surface
L_test=(M==max(M(:))); %Logic for interior voxels of surface mesh

%Construct search blurring kernel
[x,y,z]=ndgrid(-searchRadius:1:searchRadius); %kernel coordinates
d=sqrt(x.^2+y.^2+z.^2); %distance metric
k=(d<=searchRadius); %Define kernel based on radius
k=k./sum(k(:)); %Normalize kernel

%Convolve interior image with kernel
ML=convn(double(L_test),k,'same');
[~,indInternal]=max(ML(:)); %Kernel should yield max at "deep" (related to search radius) point
[I_in,J_in,K_in]=ind2sub(size(M),indInternal); %Convert to subscript coordinates

%Derive point coordinates
V_in=zeros(1,3);
[V_in(1),V_in(2),V_in(3)]=im2cart(I_in,J_in,K_in,voxelSize*ones(1,3));
V_in=V_in+imOrigin(ones(size(V_in,1),1),:);

%% Plotting image, mesh and inner point
if plotOn==1    
    figColor='w'; figColorDef='white';
    fontSize=20;
    markerSize1=round(voxelSize*25);
    faceAlpha1=0.7;
    faceAlpha2=0.2;
    
    hf=figuremax(figColor,figColorDef);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hold on;
    
    hp= patch('Faces',F,'Vertices',V,'FaceColor','g','EdgeColor','none','FaceAlpha',faceAlpha2);
    plotV(V_in,'r.','MarkerSize',markerSize1);
    
    L_plot=false(size(M));
    L_plot(:,:,K_in)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,double(M),'sk');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm-imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    L_plot=false(size(M));
    L_plot(I_in,:,:)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,M,'si');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm-imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    L_plot=false(size(M));
    L_plot(:,J_in,:)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,M,'sj');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm-imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    colormap(gray(3));
    caxis([0 2]);
    hc=colorbar;
    set(hc,'YTick',[1/3 1 5/3]);
    set(hc,'YTickLabel',{'Exterior','Boundary','Intertior'});
    
    axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
    set(gca,'FontSize',fontSize);
    drawnow;    
end

