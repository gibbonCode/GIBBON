function [N12,G12]=mergeImageData(N1,G1,N2,G2,nFac,plotOn)

% function [N12,G12]=mergeImageData(N1,G1,N2,G2,nFac,plotOn)
% ------------------------------------------------------------------------
%
% nFac=Factor for sub-division, smallest voxel dimension in v1 is aimed to be equivalent to nFac*v2_sub
%
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 2014/01/28
%
% 2014/01/28, fixed sparse input such that input sizes are consistent i.e.
% added (:) to form column
%------------------------------------------------------------------------

%%

%Creating subImage indices using first image
nSub=round(G2.v./min(G1.v(:)/nFac))'; %Sub-division specification, can be inhomogeneous since it depends on v2
M2=N2(:,:,:,1);
siz_M2_sub=size(M2).*nSub;
num_M2_sub=prod(siz_M2_sub);
[~,linearIndSub]=subImage(M2,nSub);
v2_sub=G2.v./nSub';

% if size(N1,4)==size(N2,4)
    siz4=size(N1,4);
% else
%     warning('4th dimension of images does not match, using lowest dimension');
%     siz4=min([size(N1,4) size(N2,4)]);
% end

N12=[];
disp('<<<< mergeImageData >>>>');
t_step_c=0;
for qd=1:1:siz4
    tic; %Start clock for iteration time keeping
    disp(['------ MERGING 4th DIM. STEP: ',num2str(qd),' OF ',num2str(siz4),' ------']);
    
    %Get image sets
    M1=N1(:,:,:,qd);
    M2=N2(:,:,:,qd);
        
    M1c=M1(:); %Column (helps cope with singular dimensions)
    M2c=M2(:); %Column (helps cope with singular dimensions)    
    
    %Subdeviding image set
    disp('Creating subdivided image set');
    M2_sub=M2c(linearIndSub); %Create subdevided image           
    M2_sub_c=M2_sub(:); %Column (helps cope with singular dimensions)
    
    %Coordinate grids
    disp('Defining coordinate grids');
    [I2s,J2s,K2s]=ndgrid(1:size(M2_sub,1),1:size(M2_sub,2),1:size(M2_sub,3)); %image coordinates of M2s voxel centres
    [X2s,Y2s,Z2s]=im2cart(I2s,J2s,K2s,v2_sub); %Convert above to cartesian coordinates in the colinear cartesian system of M2s
    [I2ss,J2ss,K2ss]=cart2im(X2s,Y2s,Z2s,G2.v); %Convert cartesian coordinates of M2s to image coordinates of M2
    [X2sr,Y2sr,Z2sr]=im2MRcart(I2ss,J2ss,K2ss,G2.v,G2.OR,G2.r,G2.c); %Convert to "real world" cartesian coordinates
    [I2s1,J2s1,K2s1]=MRcart2im(X2sr,Y2sr,Z2sr,G1.v,G1.OR,G1.r,G1.c); %Convert cartesian coordinates of set 2 to image coordinates of set 1
    
    disp('Mapping intensities into new image');
    %Rounding coordinates "so they snap" to voxel centre coordinates in set 1
    I2s1_r=round(I2s1); J2s1_r=round(J2s1); K2s1_r=round(K2s1);
    
    %Creat new M12 image that covers entire merged space
    max_I=max([I2s1_r(:); size(M1,1)]);
    max_J=max([J2s1_r(:); size(M1,2)]);
    max_K=max([K2s1_r(:); size(M1,3)]);
    min_I=min([I2s1_r(:); 1]);
    min_J=min([J2s1_r(:); 1]);
    min_K=min([K2s1_r(:); 1]);
    
    %Fixing potential negative indices
    I2s1_r_f=(I2s1_r-min_I)+1;
    J2s1_r_f=(J2s1_r-min_J)+1;
    K2s1_r_f=(K2s1_r-min_K)+1;
    
    %Preparing indices of M1 intensities for M12
    [I1,J1,K1]=ndgrid(1:size(M1,1),1:size(M1,2),1:size(M1,3)); %image coordinates of M1 voxel centres
    I1_f=(I1-min_I)+1;
    J1_f=(J1-min_J)+1;
    K1_f=(K1-min_K)+1;
    
    max_I_f=max([I2s1_r_f(:); I1_f(:)]);
    max_J_f=max([J2s1_r_f(:); J1_f(:)]);
    max_K_f=max([K2s1_r_f(:); K1_f(:)]);
    
    M12=zeros(max_I_f,max_J_f,max_K_f); %Initialize M12
    M12_1=zeros(max_I_f,max_J_f,max_K_f); %Initialize M12_1
    M12_2=zeros(max_I_f,max_J_f,max_K_f); %Initialize M12_2
    
    IND1_f=sub2ind(size(M12),I1_f,J1_f,K1_f);
    M12_1(IND1_f(:))=M1c;
    
    L_invalid_1=true(size(M12));
    L_invalid_1(IND1_f(:))=0;
    
    %Convert subscript to linear indices.
    IND2s1_r_f=sub2ind(size(M12),I2s1_r_f,J2s1_r_f,K2s1_r_f);
    IND2s1_r_f_uni=unique(IND2s1_r_f);
    
    % Create sparse array
    disp('Sparse array based volume averaging');
    tempSparseShift=0;%min(M2_sub_c(:))-1;
    % S = sparse(i,j,s,m,n,nzmax)
    sparseValues=M2_sub_c-tempSparseShift;
    
    intensityMap_M2_sub_M12 = sparse(1:1:numel(M2_sub),IND2s1_r_f(:),sparseValues,numel(M2_sub),numel(M12),numel(IND2s1_r_f));
        
    logicMap_M2_sub_M12 = sparse(1:1:numel(M2_sub),IND2s1_r_f(:),1,numel(M2_sub),numel(M12),numel(IND2s1_r_f));
    
    summedIntensities_M2_sub_M12 = full(sum(intensityMap_M2_sub_M12,1))';
    summedLogic_M2_sub_M12 = full(sum(logicMap_M2_sub_M12,1))'; %Number of M2_sub voxels found within a M12 voxel
    
    volumeRatio_M2_sub_M12=(summedLogic_M2_sub_M12.*prod(v2_sub))./prod(G1.v);
    % volumeRatio_M2_sub_M12=volumeRatio_M2_sub_M12./(max(volumeRatio_M2_sub_M12(:)));
    
    M12_volumeRatio_M2_sub_M12=reshape(volumeRatio_M2_sub_M12,size(M12));
    meanIntensity_M2_sub_M12=summedIntensities_M2_sub_M12./summedLogic_M2_sub_M12;
    meanIntensity_M2_sub_M12(isnan(meanIntensity_M2_sub_M12))=0;
    newIntensity_M2_sub_M12=meanIntensity_M2_sub_M12.*volumeRatio_M2_sub_M12;
    newIntensity_M2_sub_M12=reshape(newIntensity_M2_sub_M12,size(M12))+tempSparseShift;
    
    %Scale maximum if increased
    if max(newIntensity_M2_sub_M12(:))>max(M2_sub(:))
        newIntensity_M2_sub_M12=(newIntensity_M2_sub_M12./max(newIntensity_M2_sub_M12(:))).*max(M2_sub(:));
    end
    M12_2=newIntensity_M2_sub_M12;
    
    M12=(M12_1+M12_2.*M12_volumeRatio_M2_sub_M12)./(~L_invalid_1+M12_volumeRatio_M2_sub_M12);


    
    %Create new geometry parameters
    disp('Defining geometry parameters');
    v12=G1.v;
    r12=G1.r;
    c12=G1.c;
    OR12=G1.OR;
    I_OR12=1+(min_I-1);
    J_OR12=1+(min_J-1);
    K_OR12=1+(min_K-1);
    [OR12(1),OR12(2),OR12(3)]=im2MRcart(I_OR12,J_OR12,K_OR12,v12,OR12,r12,c12); %Convert image to cartesian coordinates
    
    %Collect in structure
    G12.v=v12;
    G12.r=r12;
    G12.c=c12;
    G12.OR=OR12([2 1 3]);
    
    %%
    
    if qd==1 %After first iteration
        N12=M12;
        N12=repmat(N12,[1 1 1 siz4]);
    else        
        N12(:,:,:,qd)=M12;
        if plotOn==1
            close([hf1 hf2 hf3 hf4]);
        end
    end   
    
    %%
        
    if plotOn==1
        
        %% Plot settings
        figColor='w'; figColorDef='white';
        fontSize=20;
        cMap=gray(250);
        faceAlpha1=1;
        faceAlpha2=0.5;
        edgeColor1='r';
        lineWidth1=1;
        edgeColor2='g';
        lineWidth2=1;
        
        %% PLOTTING DATA IN REAL WORLD SPACE
        
        hf1=cFigure;  set(gcf,'renderer','OpenGL');
        xlabel('X (mm)','FontSize',fontSize); ylabel('Y (mm)','FontSize',fontSize); zlabel('Z (mm)','FontSize',fontSize);
        title('Real-world cartesian space','FontSize',fontSize);
        hold on;
        
        [F,V,C]=ind2patch(true(size(M1)),M1,'vb'); %Creating patch data for selection of low voxels
        C(isnan(C))=0;
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G1.v,G1.OR,G1.r,G1.c); %Convert image to cartesian coordinates
        hp1= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        [F,V,C]=ind2patch(true(size(M2)),M2,'vb'); %Creating patch data for selection of low voxels
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G2.v,G2.OR,G2.r,G2.c); %Convert image to cartesian coordinates
        hp2= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2,'LineWidth',lineWidth2);
        
        axis equal; view(3); axis tight;  grid on;
        colormap(cMap); colorbar;
        set(gca,'FontSize',fontSize);
        drawnow;
        
        %% PLOTTING DATA IN IMAGE SPACE OF SET 1
        
        hf2=cFigure; set(gcf,'renderer','OpenGL');
        xlabel('J','FontSize',fontSize); ylabel('I','FontSize',fontSize); zlabel('K','FontSize',fontSize);
        title('Image space of set 1','FontSize',fontSize);
        hold on;
        
        [F,V,C]=ind2patch(true(size(M1)),M1,'vb'); %Creating patch data for selection of low voxels
        C(isnan(C))=0;
        hp1= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        [F,V,C]=ind2patch(true(size(M2)),M2,'vb'); %Creating patch data for selection of low voxels
        C(isnan(C))=0;
        [Xs,Ys,Zs]=im2MRcart(V(:,2),V(:,1),V(:,3),G2.v,G2.OR,G2.r,G2.c); %Convert image to cartesian coordinates
        [V(:,2),V(:,1),V(:,3)]=MRcart2im(Xs,Ys,Zs,G1.v,G1.OR,G1.r,G1.c); %Convert cartesian coordinates to image coordinates of set 1
        hp2= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor2,'FaceAlpha',faceAlpha2,'LineWidth',lineWidth2);
        
        axis equal; view(3); axis tight;  grid on;
        colormap(cMap); colorbar;
        set(gca,'FontSize',fontSize);
        drawnow;
        
        %% PLOTTING MERGED IMAGE IN IMAGE SPACE
        
        hf3=cFigure; set(gcf,'renderer','OpenGL');
        xlabel('J','FontSize',fontSize); ylabel('I','FontSize',fontSize); zlabel('K','FontSize',fontSize);
        title('Image space of merged set','FontSize',fontSize);
        hold on;
        
        [F,V,C]=ind2patch(~isnan(M12),M12,'vb'); %Creating patch data for selection of low voxels
        C(isnan(C))=0;
        hp= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        axis equal; view(3); axis tight;  grid on;
        colormap(cMap); colorbar; %caxis([0 1]);
        set(gca,'FontSize',fontSize);
        drawnow;

        %%
        
        hf4=cFigure;  set(gcf,'renderer','OpenGL');
        xlabel('J','FontSize',fontSize); ylabel('I','FontSize',fontSize); zlabel('K','FontSize',fontSize);
        title('Image space of set 1','FontSize',fontSize);
        hold on;
        
        L_slice=false(size(M12)); L_slice(round(size(M12,1)/2),:,:)=1; L_slice(isnan(M12))=0;
        [F,V,C]=ind2patch(L_slice,M12,'si'); %Creating patch data for selection of low voxels
%         C(isnan(C))=0;
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G12.v,G12.OR,G12.r,G12.c); %Convert image to cartesian coordinates
        hp= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        L_slice=false(size(M12)); L_slice(:,round(size(M12,2)/2),:)=1; L_slice(isnan(M12))=0;
        [F,V,C]=ind2patch(L_slice,M12,'sj'); %Creating patch data for selection of low voxels
%         C(isnan(C))=0;
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G12.v,G12.OR,G12.r,G12.c); %Convert image to cartesian coordinates
        hp= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        L_slice=false(size(M12)); L_slice(:,:,round(size(M12,3)/2))=1; L_slice(isnan(M12))=0;
        [F,V,C]=ind2patch(L_slice,M12,'sk'); %Creating patch data for selection of low voxels
%         C(isnan(C))=0;
        [V(:,1),V(:,2),V(:,3)]=im2MRcart(V(:,2),V(:,1),V(:,3),G12.v,G12.OR,G12.r,G12.c); %Convert image to cartesian coordinates
        hp= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor',edgeColor1,'FaceAlpha',faceAlpha1,'LineWidth',lineWidth1);
        
        axis equal; view(3); axis tight;  grid on;
        colormap(cMap); colorbar; %caxis([0 1]);
        set(gca,'FontSize',fontSize);
        drawnow;                
                
    end
    
    %% Time keeping
    t_step=toc;
    t_step_c=t_step_c+t_step;
    t_left_min=(t_step.*(siz4-qd))./60;
    disp('Finished current dynamic');
    disp(['Time taken: ',num2str(t_step),' s, Estimated time left: ',num2str(t_left_min),' min.']);    
    
end

% if plotOn==1
%     close([hf1 hf2 hf3 hf4]);
% end
 
%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2019  Kevin Mattheus Moerman
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
