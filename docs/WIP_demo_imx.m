%% DEMO_imx
% Below is a demonstration of the features of the |triSurf2Im| function

%%
clear; close all; clc;

%%
% Plot settings
fontSize=10;
faceAlpha1=1;
faceAlpha2=0.3;

%% Simulate image
% Defining an example triangulated surface model

% Defining a deformed and rotated torus shape
r=1; %Sphere radius
rc=2; %Central radius
nr=16;
nc=30;
ptype='tri';
[F,V]=patchTorus(r,nr,rc,nc,ptype);
[THETA,RHO] = cart2pol(V(:,1),V(:,2));

%%
% Setting control parameters

% Defining the full set of possible control parameters
voxelSize=0.15; % The output image voxel size.
imOrigin=min(V,[],1)-2*voxelSize;
imMax=max(V,[],1)+2*voxelSize;
imSiz=round((imMax-imOrigin)/voxelSize);
imSiz=imSiz([2 1 3]); %Image size (x, y corresponds to j,i in image coordinates, hence the permutation)

% Using |triSurf2Im| function to convert patch data to image data
[M,G,bwLabels]=triSurf2Im(F,V,voxelSize,imOrigin,imSiz);

%%
% Plotting the results

hf1=cFigure;
subplot(1,2,1);
title('Closed triangulated surface','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;
patch('Faces',F,'Vertices',V,'FaceColor','g','EdgeColor','k','FaceAlpha',faceAlpha1);
camlight('headlight'); lighting flat;
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);

subplot(1,2,2);
title('Boundary, intertior and exterior image','FontSize',fontSize);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

patch('Faces',F,'Vertices',V,'FaceColor','g','EdgeColor','none','FaceAlpha',faceAlpha2);
L_plot=false(size(M));
L_plot(:,:,round(size(M,3)/2))=1;
[Fm,Vm,Cm]=ind2patch(L_plot,double(M),'sk');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);

L_plot=false(size(M));L_plot(round(size(M,1)/2),:,:)=1;
[Fm,Vm,Cm]=ind2patch(L_plot,M,'si');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);

L_plot=false(size(M));L_plot(:,round(size(M,2)/2),:)=1;
[Fm,Vm,Cm]=ind2patch(L_plot,M,'sj');
[Vm(:,1),Vm(:,2),Vm(:,3)]=im2cart(Vm(:,2),Vm(:,1),Vm(:,3),voxelSize*ones(1,3));
Vm=Vm+imOrigin(ones(size(Vm,1),1),:);
patch('Faces',Fm,'Vertices',Vm,'FaceColor','flat','CData',Cm,'EdgeColor','k','FaceAlpha',faceAlpha1);

colormap(gray(3)); caxis([0 2]);
hc=colorbar;
set(hc,'YTick',[1/3 1 5/3]);
set(hc,'YTickLabel',{'Exterior','Boundary','Intertior'});
axis equal; view(3); axis tight;  grid on;  set(gca,'FontSize',fontSize);
drawnow;

%%
% Start segmentation using |imx|

v=voxelSize*ones(1,3);
hf=imx(M,v);
drawnow; hold on;

Vt=V-imOrigin(ones(size(V,1),1),:);
gpatch(F,Vt,'g','none',0.5);

% msgbox('Try to load in imseg_torus_1.mat or imseg_torus_2.mat from...GIBBON/data/temp/imseg/','Hint','help')

%%
% Turn contours into segmentation

pathName='/mnt/data/MATLAB/GIT/GIBBON/data/temp/imseg/';
loadNames={'imseg_torus_1','imseg_torus_2'};

close all;

hf1q=sv3(M,v);
drawnow;

numContours=numel(loadNames);
numSlices=size(M,3);

plotColors=gjet(4);

levelSetType=2;
vizStruct.colormap=gjet(250);

%%

KK=repmat(M,1,1,1,2);
for q=1:1:2
    loadName=fullfile(pathName,loadNames{q});
    load(loadName);
    Vcs=saveStruct.ContourSet;
    [K]=contour2levelset(M,v,Vcs,levelSetType);
    KK(:,:,:,q)=K;
    
    for qSlice=1:1:numSlices
        numSubContours=numel(Vcs{qSlice});
        for qSub=1:1:numSubContours
            Vd=Vcs{qSlice}{qSub}; %Current contour
            if ~isempty(Vd)
                hp=plotV(Vd,'w-');
                hp.Color=plotColors(q,:);
                hp.LineWidth=5;
            end
        end
    end    
end
gpatch(F,Vt,'r','none',0.5);
drawnow;

%%

KK(:,:,:,2)=-KK(:,:,:,2);
K=max(KK,[],4);%KK(:,:,:,1);
K(isnan(K))=nanmax(K(:));

hf1q=sv3(K,v,vizStruct);
% 
% hf2q=sv3(K2,v,vizStruct);

%%

controlPar.contourLevel=0;
controlPar.voxelSize=v;
controlPar.nSub=[1 1 1];
controlPar.capOpt=1;

[Fi,Vi]=levelset2isosurface(K,controlPar);
Fi=fliplr(Fi); %Invert orientation
[Fi,Vi]=triSurfRemoveThreeConnect(Fi,Vi,[]);

cParSmooth.Method='HC';
cParSmooth.Alpha=0.1;
cParSmooth.Beta=0.5;
cParSmooth.n=15;
[Vi]=patchSmooth(Fi,Vi,[],cParSmooth);

%%

% Di=minDist(Vi,Vt);
[Di]=triSurfSetDist(Fi,Vi,F,Vt,'dist-ray');

%%
cFigure; 
hold on; 
gpatch(F,Vt,0.5*ones(1,3),'none',0.5);
hp=gpatch(Fi,Vi,Di,'none',1);
caxis([min(Di(:)) max(Di(:))]);
colormap(gjet(250)); colorbar; 
axisGeom;
drawnow

%%





