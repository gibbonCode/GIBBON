function [V_inner]=getInnerPoint(F_input,V_input,searchRadius,voxelSize,plotOn)

%%

if isa(F_input,'cell')
    [F_input,V_input,C_input]=joinElementSets(F_input,V_input);
else
    C_input=ones(size(F_input,1),1);
end

%% Get an inner point

%Convert patch data to image
[M,G,~]=patch2Im(F_input,V_input,C_input,voxelSize);
L_test=(M==1); %Logic for interior voxels of surface mesh
M=double(L_test); 
imOrigin=G.origin; %Image origin can be used to allign image with surface


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
V_inner=zeros(1,3);
[V_inner(1),V_inner(2),V_inner(3)]=im2cart(I_in,J_in,K_in,voxelSize*ones(1,3));
V_inner=V_inner+imOrigin(ones(size(V_inner,1),1),:);

%% Plotting image, mesh and inner point
if plotOn==1    
    figColor='w'; figColorDef='white';
    fontSize=20;
    markerSize1=50;%round(voxelSize*25);
    faceAlpha1=0.7;
    faceAlpha2=0.2;
    
    hf=figuremax(figColor,figColorDef);
    xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
    hold on;
    
    hp= patch('Faces',F_input,'Vertices',V_input,'FaceColor','flat','CData',C_input,'EdgeColor','none','FaceAlpha',faceAlpha2);
    plotV(V_inner,'k.','MarkerSize',markerSize1);
    
    L_plot=false(size(M));
    L_plot(:,:,K_in)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,double(M),'sk');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor',0.5*ones(1,3),'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    L_plot=false(size(M));
    L_plot(I_in,:,:)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,M,'si');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor',0.5*ones(1,3),'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    L_plot=false(size(M));
    L_plot(:,J_in,:)=1;
    L_plot=L_plot&M>0;
    [Fm,Vm,Cm]=ind2patch(L_plot,M,'sj');
    [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
    Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
    patch('Faces',Fm,'Vertices',Vm,'FaceColor',0.5*ones(1,3),'EdgeColor','k','FaceAlpha',faceAlpha1);
    
    colormap(gjet(numel(unique(C_input(:))))); colorbar; 
    
    axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
    set(gca,'FontSize',fontSize);
    drawnow;    
end

