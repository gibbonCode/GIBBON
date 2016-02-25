function hf=sv3(varargin)

%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
    case 2
        M=varargin{1};
        v=varargin{2};  
end

M=double(M);

%%
% Plot settings
fontColor='w';
fontSize=20;
cMap=gray(250);

figStruct.Name='GIBBON: Slice viewer'; %Figure name
figStruct.Color='k'; %Figure background color
figStruct.ColorDef='black'; %Setting colordefinitions to black
% figStruct.ScreenOffset=20; %Setting spacing of figure with respect to screen edges

%%

%Defining row, column and slice indicices for slice patching
sliceIndexI=round(size(M,1)/2); %(close to) middle row
sliceIndexJ=round(size(M,2)/2); %(close to) middle column
sliceIndexK=round(size(M,3)/2); %(close to) middle slice

Tf=[0 100]; %Threshold

nTickMajor=20;
tickSizeMajor_I=round(size(M,1)/nTickMajor);
tickSizeMajor_J=round(size(M,2)/nTickMajor);
tickSizeMajor_K=round(size(M,3)/nTickMajor);

%%

[ax,ay,az]=im2cart([size(M,1)+1 0],[size(M,2)+1 0],[size(M,3)+1 0],v);

hf=cFigure(figStruct);
axis equal; axis tight; view(3);  axis vis3d; axis([ax(2) ax(1) ay(2) ay(1) az(2) az(1)]); grid on; box on; hold on;
colormap(cMap); colorbar;
caxis([min(M(:)) max(M(:))]);
set(gca,'fontSize',fontSize);
drawnow;

w=50; %Scrollbar width
jSlider_I = javax.swing.JSlider(1,size(M,1));
javacomponent(jSlider_I,[0,0,w,round(hf.Position(4))]);
set(jSlider_I, 'MajorTickSpacing',tickSizeMajor_I, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_I,1}},'Orientation',jSlider_I.VERTICAL);

jSlider_J = javax.swing.JSlider(1,size(M,2));
javacomponent(jSlider_J,[1*w,0,w,round(hf.Position(4))]);
set(jSlider_J, 'MajorTickSpacing',tickSizeMajor_J, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_J,2}},'Orientation',jSlider_J.VERTICAL);

jSlider_K = javax.swing.JSlider(1,size(M,3));
javacomponent(jSlider_K,[2*w,0,w,round(hf.Position(4))]);
set(jSlider_K, 'MajorTickSpacing',tickSizeMajor_K, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,jSlider_K,3}},'Orientation',jSlider_K.VERTICAL);

jSlider_T = com.jidesoft.swing.RangeSlider(0,100,Tf(1),Tf(2));  % min,max,low,high
javacomponent(jSlider_T,[3*w,0,w,round(hf.Position(4))]);
set(jSlider_T, 'MajorTickSpacing',25, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
    'Background',java.awt.Color.white, 'snapToTicks',false, 'StateChangedCallback',{@setThreshold,{hf,jSlider_T,jSlider_I,jSlider_J,jSlider_K}},'Orientation',jSlider_T.VERTICAL);

%% Set resize function 

set(hf,'ResizeFcn',{@setScrollSizeFunc,{hf,w,jSlider_T,jSlider_I,jSlider_J,jSlider_K}});

%%
hf.UserData.M=M;
hf.UserData.v=v;
hf.UserData.patchTypes={'si','sj','sk'};
hf.UserData.sliceIndices=[sliceIndexI sliceIndexJ sliceIndexK];
hf.UserData.fontColor=fontColor;
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
hf.UserData.sliceIndices(dirOpt)=sliceIndex;
sliceIndices=hf.UserData.sliceIndices;

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

title(['IJK-view ',num2str(sliceIndices(1)),' ',num2str(sliceIndices(2)),' ',num2str(sliceIndices(3))],'color',hf.UserData.fontColor);

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

function setScrollSizeFunc(~,~,inputCell)
hf=inputCell{1};
w=inputCell{2};

for q=3:6
    jSlider=inputCell{q};
    javacomponent(jSlider,[w*(q-3),0,w,round(hf.Position(4))]);
end

end

