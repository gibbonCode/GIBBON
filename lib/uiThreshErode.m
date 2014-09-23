function [L_BG]=uiThreshErode(M,thresholdInitial,blurKernelSize,groupCropOption)

%% PLOT SETTINGS
figColor='w'; fig_colordef='white';
fontSize=20;
cMap=autumn(250);
falpha=1;
patchTypes={'sj','si','sk','v'};
ptype=3;
sliceScale=2;
noSlices=7;
S=round(linspace(1,size(M,3),noSlices)); %Slice indices for plotting

%%
M_original=M;
M=double(M);

%% Normalizing the image data
M=dataNorm(M,1);

%% APPLY PRE-BLURRING
if blurKernelSize>0
    %Defining erosion/dilation kernel
    k=blurKernelSize;
    p=k-round(k./2);
    
    hb=gauss_kernel(k,ndims(M),1.5,'width'); %Gaussian kernel
       
    M_rep=zeros(size(M)+(2.*p));
    M_rep(p+(1:size(M,1)),p+(1:size(M,2)),p+(1:size(M,3)))=M;
    M = convn(M_rep,hb,'valid');
end

%% Defining erosion/dilation kernel
k=3;
hg=ones(k,k,k); 
p=k-round(k./2);

%% Start loop
done=0;
while done==0
    %Get initial threshold level    
    T_threshold=thresholdInitial; %default start value
    done=0;
    qc=1;
    runMode=1;
    while done==0
        
        switch runMode
            case 1 %Determining threshold
                L_BG=M>T_threshold;
            case 2 %Setting erosion/dilation
                %"Blurring current mask"
                L_BG_rep=zeros(size(M)+(2.*p));
                L_BG_rep(p+(1:size(M,1)),p+(1:size(M,2)),p+(1:size(M,3)))=L_BG;
                L_BG_blur = convn(double(L_BG_rep),hg,'valid');
        end
        
        if qc==1; %Plot figure on first iteration
            
            %Cropped image
            L_slices=false(size(M));
            L_slices(:,:,S)=L_BG(:,:,S);
            IND_slices=find(L_slices); clear L_slices;
            [Fs,Vs,C_slice]=ind2patch(IND_slices,M_original,patchTypes{ptype});
            Vs(:,3)=sliceScale.*Vs(:,3);
            
            hf1=figuremax(figColor,fig_colordef);
            title(['Threshold is ',num2str(T_threshold),'*mean, press up to increase or down to decrease (by 10%), press space to keep and continue'],'FontSize',20);
            hold on; xlabel('X-J','FontSize',20);ylabel('Y-I','FontSize',20);zlabel('Z-K','FontSize',20);
            hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',C_slice,'FaceColor','flat');
            set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); colorbar; 
            set(gca,'FontSize',fontSize);
            drawnow;
            
        else %redefine patch data
            delete(hs); %remove patch data from figure
            
            %redefine patch data for the cropped image
            L_slices=false(size(M));
            L_slices(:,:,S)=L_BG(:,:,S);
            IND_slices=find(L_slices); clear L_slices;
            [Fs,Vs,C_slice]=ind2patch(IND_slices,M_original,patchTypes{ptype});
            Vs(:,3)=sliceScale.*Vs(:,3);
            
            %patch again
            if runMode==1
                title(['Threshold is ',num2str(T_threshold),'*mean, press up to increase or down to decrease (by 10%), press space to keep and continue'],'FontSize',20);
            else
                title('Current cropping press up/down to dilate/erode, WARNING EROSION MAY REMOVE ENTIRE SLICES','FontSize',20);
            end
            
            hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',C_slice,'FaceColor','flat');
            set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); colorbar; 
            set(gca,'FontSize',fontSize);
            drawnow;
        end
        
        [~,~,b]=qginput(1);
        switch b
            case 30 %  up arrow key => increase threshold by 10% OR dilate
                if runMode==1
                    T_threshold=T_threshold.*1.1;
                else
                    L_BG=(L_BG_blur>0); %blurred boundary elements that increased due to blur are now set to 1 as well so added to mask
                end
            case 31 %down arrow key => decrease threshold by 10% OR erode
                if runMode==1
                    T_threshold=T_threshold.*0.9;
                else
                    L_BG=(L_BG_blur==sum(hg(:))); %blurred boundary elements that decreased due to blur are now set to 0 as well so removed from mask
                end
            case 32 %spacebar key => confirms current threshold result, exits while loop
                if runMode==1
                    runMode=2;
                else
                    done=1;
                end
        end
        qc=qc+1;
    end
    
    %% GROUPING BASED CROPPING
    if groupCropOption==1
        
        GROUP_STRUCT = bwconncomp(L_BG,6);
        IMAGE_OBJECTS=GROUP_STRUCT.PixelIdxList;
        
        %Filtering out groups that are too small
        group_sizes = cell2mat(cellfun(@(x) max(size(x)), IMAGE_OBJECTS,'UniformOutput',0)'); %get group sizes
        [~,ind_max]=max(group_sizes); %Index of largest group
        IND_OBJECT=IMAGE_OBJECTS{ind_max}; %indices of voxels for largest group
        
        %Redefining background logic using group
        L_BG=false(size(L_BG)); 
        L_BG(IND_OBJECT)=1;   
        
        %Plotting the result
        delete(hs); %remove patch data from figure
        
        %redefine patch data for the cropped image
        L_slices=false(size(M));
        L_slices(:,:,S)=L_BG(:,:,S);
        IND_slices=find(L_slices); clear L_slices;
        [Fs,Vs,C_slice]=ind2patch(IND_slices,M_original,patchTypes{ptype});
        Vs(:,3)=sliceScale.*Vs(:,3);
                
        title('Grouping result','FontSize',20);
               
        hs=patch('Faces',Fs,'Vertices',Vs,'EdgeColor','none', 'CData',C_slice,'FaceColor','flat');
        set(hs,'FaceAlpha',falpha); hold on; grid on;  view(3); axis equal; axis tight; axis vis3d;  colormap(cMap); 
        set(gca,'FontSize',fontSize);
        drawnow;
    end

    %%
    Q = questdlg('Would you like to keep this result?','Keep result?','YES','NO','YES');
    switch Q
        case 'YES'
            done=1;
        case 'NO'
            done=0;
    end
    close(hf1)
    
end
