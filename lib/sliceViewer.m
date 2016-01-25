function hf=sliceViewer(M,G)

%% Parse input

v=G.v;
M=double(M);

%%
% Plot settings
fontSize=20;
cMap=gray(250);

figStruct.Name='GIBBON: Slice viewer'; %Figure name
figStruct.Color='w'; %Figure background color
figStruct.ColorDef='white'; %Setting colordefinitions to black
% figStruct.ScreenOffset=20; %Setting spacing of figure with respect to screen edges

%%

%Defining row, column and slice indicices for slice patching
sliceIndexI=round(size(M,1)/2); %(close to) middle row
sliceIndexJ=round(size(M,2)/2); %(close to) middle column
sliceIndexK=round(size(M,3)/2); %(close to) middle slice

Tf=[0 100]; %Threshold

%%

[ax,ay,az]=im2cart([size(M,1)+1 0],[size(M,2)+1 0],[size(M,3)+1 0],v);

%Open figure
hf=cFigure(figStruct);
% title('MRI visualisation, slices and voxels in cartesian coordinates with aid of voxel size');

axis equal; axis tight; view(3);  axis vis3d; axis([ax(2) ax(1) ay(2) ay(1) az(2) az(1)]); grid on; box on;
colormap(cMap); colorbar;
caxis([min(M(:)) max(M(:))]);
% colorbar('Ticks',[linspace(min(M(:)),max(M(:)),10)]);
% caxis manual;
set(gca,'fontSize',fontSize);
drawnow;

w=50;
jSlider_I = javax.swing.JSlider(1,size(M,1));% jRangeSlider = javacomponent(jRangeSlider, [0,0,200,80], gcf);
javacomponent(jSlider_I,[0,0,w,round(hf.Position(4)/2)]);
set(jSlider_I, 'MajorTickSpacing',4, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_I,1}},'Orientation',jSlider_I.VERTICAL);
% set(jSlider_I,'Opaque',false);
set(jSlider_I,'Foreground',java.awt.Color.black);


jSlider_J = javax.swing.JSlider(1,size(M,2));% jRangeSlider = javacomponent(jRangeSlider, [0,0,200,80], gcf);
javacomponent(jSlider_J,[1*w,0,w,round(hf.Position(4)/2)]);
set(jSlider_J, 'MajorTickSpacing',4, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_J,2}},'Orientation',jSlider_I.VERTICAL);

jSlider_K = javax.swing.JSlider(1,size(M,3));% jRangeSlider = javacomponent(jRangeSlider, [0,0,200,80], gcf);
javacomponent(jSlider_K,[2*w,0,w,round(hf.Position(4)/2)]);
set(jSlider_K, 'MajorTickSpacing',4, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_K,3}},'Orientation',jSlider_I.VERTICAL);

% jSlider_T = javax.swing.JSlider(0,100);% jRangeSlider = javacomponent(jRangeSlider, [0,0,200,80], gcf);
% javacomponent(jSlider_T,[3*w,0,w,round(hf.Position(4)/2)]);
% set(jSlider_T, 'MajorTickSpacing',5, 'MinorTickSpacing',25, 'PaintTicks',true, 'PaintLabels',true,...
%     'Background',java.awt.Color.white, 'snapToTicks',false, 'StateChangedCallback',{@setThreshold,{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K}},'Orientation',jSlider_I.VERTICAL);

jSlider_T = com.jidesoft.swing.RangeSlider(0,100,Tf(1),Tf(2));  % min,max,low,high
javacomponent(jSlider_T,[3*w,0,w,round(hf.Position(4)/2)]);
set(jSlider_T, 'MajorTickSpacing',25, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',false, 'StateChangedCallback',{@setThreshold,{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K}},'Orientation',jSlider_I.VERTICAL);

%%
hf.UserData.M=M;
hf.UserData.v=v;
hf.UserData.patchTypes={'si','sj','sk'};

%%
set(jSlider_I,'Value',sliceIndexI);
set(jSlider_J,'Value',sliceIndexJ);
set(jSlider_K,'Value',sliceIndexK);
set(jSlider_T,'HighValue',Tf(2));
set(jSlider_T,'LowValue',Tf(1));

setThreshold([],[],{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K});

end

function plotSlice(~,~,inputCell)

hf=inputCell{1};
jSlider=inputCell{2};
dirOpt=inputCell{3};
sliceIndex = get(jSlider,'Value');

M=hf.UserData.M;
logicThreshold=hf.UserData.logicThreshold;
v=hf.UserData.v;
patchType=hf.UserData.patchTypes{dirOpt}; 

logicPatch=false(size(M));
switch dirOpt
    case 1
        logicPatch(sliceIndex,:,:)=1;
    case 2
        logicPatch(:,sliceIndex,:)=1;
    case 3
        logicPatch(:,:,sliceIndex)=1;
end
logicPatch=logicPatch & logicThreshold;

[F,V,C]=ind2patch(logicPatch,M,patchType);
[V(:,1),V(:,2),V(:,3)]=im2cart(V(:,2),V(:,1),V(:,3),v);

if isfield(hf.UserData,'hp');
    try
        delete(hf.UserData.hp(dirOpt));
    catch
    end
end

hf.UserData.hp(dirOpt)= patch('Faces',F,'Vertices',V,'FaceColor','flat','CData',C,'EdgeColor','none');
drawnow;

end

function setThreshold(~,~,inputCell)

hf=inputCell{1};
jSlider_T=inputCell{2};
jSlider_I=inputCell{3};
jSlider_J=inputCell{4};
jSlider_K=inputCell{5};

Tf(1) = get(jSlider_T,'LowValue');
Tf(2) = get(jSlider_T,'HighValue');

M=hf.UserData.M;
W=max(M(:))-min(M(:));

T_low=min(M(:))+(W*Tf(1)/100);
T_high=min(M(:))+(W*Tf(2)/100);
logicThreshold=(M>=T_low & M<=T_high);
hf.UserData.logicThreshold=logicThreshold;

plotSlice([],[],{hf,jSlider_I,1});
plotSlice([],[],{hf,jSlider_J,2});
plotSlice([],[],{hf,jSlider_K,3});

end
