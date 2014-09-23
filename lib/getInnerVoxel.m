function [indInternal]=getInnerVoxel(L,searchRadius,plotOn)

if any(L(:))
    %% Get an inner voxel
    
    %Construct search blurring kernel
    [x,y,z]=ndgrid(-searchRadius:1:searchRadius); %kernel coordinates
    d=sqrt(x.^2+y.^2+z.^2); %distance metric
    k=(d<=searchRadius); %Define kernel based on radius
    k=k./sum(k(:)); %Normalize kernel
    
    %Convolve interior image with kernel
    ML=convn(double(L),k,'same');
    ML(~L)=0;
    [~,indInternal]=max(ML(:)); %Kernel should yield max at "deep" (related to search radius) point
    
    %% Plotting image, mesh and inner point
    if plotOn==1
        voxelSize=[1 1 1];
        [I_in,J_in,K_in]=ind2sub(size(L),indInternal); %Convert to subscript coordinates
        [X_in,Y_in,Z_in]=im2cart(I_in,J_in,K_in,voxelSize);
        
        %Plot settings
        figColor='w'; figColorDef='white';
        fontSize=20;
        markerSize1=round(max(voxelSize)*25);
        faceAlpha1=1;
        
        hf=figuremax(figColor,figColorDef);
        xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
        hold on;
        
        %Found voxel location
        plot3(X_in,Y_in,Z_in,'r.','MarkerSize',markerSize1);
        
        %Local slices
        L_plot=false(size(L));
        L_plot(:,:,K_in)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=ind2patch(L_plot,double(L),'sk');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
        
        L_plot=false(size(L));
        L_plot(I_in,:,:)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=ind2patch(L_plot,L,'si');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
        
        L_plot=false(size(L));
        L_plot(:,J_in,:)=1;
        L_plot=L_plot&L>0;
        [Fm,Vm,Cm]=ind2patch(L_plot,L,'sj');
        [Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize);
        patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);
        
        axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
        set(gca,'FontSize',fontSize);
        colormap(gray(3));
        drawnow;
    end
else
    indInternal=[];
end
