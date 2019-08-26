function [varargout]=sv3(varargin)

% function [hf]=sv3(M,v,optionStruct)
% ------------------------------------------------------------------------
% sv3 (slice view 3D) is a 3D slice viewer function. The 3D image data M is
% rendered using 3 mutally orthogonal slices. The option input v is the
% voxel size which scales the image appropriately in the 3 directions. The
% 3rd optional input vixOptStruct which can be used to make custom
% visualization settings. The default structure containts the following: 
%
% optionStructDefault.colormap=gray(250); %colormap
% optionStructDefault.clim=[min(M(~isnan(M))) max(M(~isnan(M)))]; %color limits
% optionStructDefault.fontColor='w'; %font color
% optionStructDefault.fontSize=20; %font size
% optionStructDefault.figStruct=figStruct; %figure options (see cFigure)
% optionStructDefault.sliceIndices=round(size(M)/2); %Default mid-slices
% optionStructDefault.updateFrequency=100; %Max update frequency during slider
% drag
%
% See also: sliceViewer, sv2, imx
%
% Change log: 
% 2018/06/06 Added initial slice indices as option to input structure
% 2018/06/06 Added basic description at the top of this function
% 2019/08/09 Changed to use uicontrol slider rather than java slider due
% to future removal of javacomponent
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1
        M=varargin{1};
        v=ones(1,3);
        optionStruct=[];
    case 2
        M=varargin{1};
        v=varargin{2}; 
        optionStruct=[];
    case 3
        M=varargin{1};
        v=varargin{2};
        optionStruct=varargin{3};
end
M=double(M); %Conver the image to a double

figStruct.Name='GIBBON: Slice viewer'; %Figure name
figStruct.Color='k'; %Figure background color
figStruct.ColorDef='black'; %Setting colordefinitions to black

optionStructDefault.colormap=gray(250); %colormap
optionStructDefault.clim=[min(M(~isnan(M))) max(M(~isnan(M)))]; %color limits
optionStructDefault.fontColor='w'; %font color
optionStructDefault.fontSize=20; %font size
optionStructDefault.figStruct=figStruct; %figure options (see cFigure)
optionStructDefault.sliceIndices=round(size(M)/2); %Default mid-slices
optionStructDefault.alphaLevel=1;
optionStructDefault.origin=[0 0 0];
optionStructDefault.updateFrequency=10;

[optionStruct]=structComplete(optionStruct,optionStructDefault,1);

M=double(M);

%%
% Plot settings
scrollBarWidth=20; %Scrollbar width
fontColor=optionStruct.fontColor;
fontSize=optionStruct.fontSize;
cMap=optionStruct.colormap;
cLim=optionStruct.clim;
figStruct=optionStruct.figStruct;
alphaLevel=optionStruct.alphaLevel;
originLoc=optionStruct.origin;
updateFrequency=optionStruct.updateFrequency; 

if diff(cLim)<eps
   cLim=cLim+[-1 1];
end

%%

%Defining row, column and slice indicices for slice patching
sliceIndexI=optionStruct.sliceIndices(1); %(close to) middle row
sliceIndexJ=optionStruct.sliceIndices(2); %(close to) middle column
sliceIndexK=optionStruct.sliceIndices(3); %(close to) middle slice

%%

[ax,ay,az]=im2cart([size(M,1)+1 0],[size(M,2)+1 0],[size(M,3)+1 0],v);

navString=['I: ',num2str(sliceIndexI),', J:  ',num2str(sliceIndexJ),', K: ',num2str(sliceIndexK)];

hf=cFigure(figStruct);
ht=gtitle(navString);
axis equal; axis tight; view(3);  axis vis3d; axis([ax(2) ax(1) ay(2) ay(1) az(2) az(1)]); grid on; box on; hold on;
colormap(cMap); colorbar;
caxis(cLim);
set(gca,'fontSize',fontSize);
drawnow;

%%

%Initialize sliders
hSlider_I= uicontrol(hf,'Style','slider','Position',[0,0,scrollBarWidth,round(hf.Position(4))]);
set(hSlider_I,'Value',sliceIndexI,'Min',1,'Max',size(M,1),'SliderStep',[1/(size(M,1)-1) 1/(size(M,1)-1)]);
hSlider_I.Callback={@plotSlice,{hf,hSlider_I,1}};
addlistener(hSlider_I,'ContinuousValueChange',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_I,1}));
addlistener(hSlider_I,'Value','PostSet',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_I,1}));

hSlider_J= uicontrol(hf,'Style','slider','Position',[1*scrollBarWidth,0,scrollBarWidth,round(hf.Position(4))]);
set(hSlider_J,'Value',sliceIndexJ,'Min',1,'Max',size(M,2),'SliderStep',[1/(size(M,2)-1) 1/(size(M,2)-1)]);
hSlider_J.Callback={@plotSlice,{hf,hSlider_J,2}};
addlistener(hSlider_J,'ContinuousValueChange',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_J,2}));
addlistener(hSlider_J,'Value','PostSet',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_J,2}));

hSlider_K= uicontrol(hf,'Style','slider','Position',[2*scrollBarWidth,0,scrollBarWidth,round(hf.Position(4))]);
set(hSlider_K,'Value',sliceIndexK,'Min',1,'Max',size(M,3),'SliderStep',[1/(size(M,3)-1) 1/(size(M,3)-1)]);
hSlider_K.Callback={@plotSlice,{hf,hSlider_K,3}};
addlistener(hSlider_K,'ContinuousValueChange',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_K,3}));
addlistener(hSlider_K,'Value','PostSet',@(hObject, event) plotSlice(hObject,event,{hf,hSlider_K,3}));

%%
 
% hSlider_I = javax.swing.JSlider(1,size(M,1));
% javacomponent(hSlider_I,[0,0,scrollBarWidth,round(hf.Position(4))]);
% set(hSlider_I, 'MajorTickSpacing',tickSizeMajor_I, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
%     'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,hSlider_I,1}},'Orientation',hSlider_I.VERTICAL);
% 
% hSlider_J = javax.swing.JSlider(1,size(M,2));
% javacomponent(hSlider_J,[1*scrollBarWidth,0,scrollBarWidth,round(hf.Position(4))]);
% set(hSlider_J, 'MajorTickSpacing',tickSizeMajor_J, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
%     'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,hSlider_J,2}},'Orientation',hSlider_J.VERTICAL);
% 
% hSlider_K = javax.swing.JSlider(1,size(M,3));
% javacomponent(hSlider_K,[2*scrollBarWidth,0,scrollBarWidth,round(hf.Position(4))]);
% set(hSlider_K, 'MajorTickSpacing',tickSizeMajor_K, 'MinorTickSpacing',1, 'PaintTicks',true, 'PaintLabels',true,...
%     'Background',java.awt.Color.white, 'snapToTicks',true, 'StateChangedCallback',{@plotSlice,{hf,hSlider_K,3}},'Orientation',hSlider_K.VERTICAL);

%% Set resize function 

set(hf,'ResizeFcn',{@setScrollSizeFunc,{hf,scrollBarWidth,hSlider_I,hSlider_J,hSlider_K}});

%%

t=clock;
hf.UserData.sv3.time=t;
hf.UserData.sv3.Name=figStruct.Name;
hf.UserData.sv3.M=M;
hf.UserData.sv3.v=v;
hf.UserData.sv3.patchTypes={'si','sj','sk'};
hf.UserData.sv3.sliceIndices=[sliceIndexI sliceIndexJ sliceIndexK];
hf.UserData.sv3.fontColor=fontColor;
hf.UserData.sv3.fontSize=fontSize;
hf.UserData.sv3.hp=nan(3,1);
hf.UserData.sv3.ht=ht;
hf.UserData.sv3.M_plot=M;
hf.UserData.sv3.alphaLevel=alphaLevel;
hf.UserData.sv3.origin=originLoc;
hf.UserData.sv3.sliderHandles=[hSlider_I,hSlider_J,hSlider_K];
hf.UserData.sv3.updateFrequency=updateFrequency;

%%
set(hSlider_I,'Value',sliceIndexI);
set(hSlider_J,'Value',sliceIndexJ);
set(hSlider_K,'Value',sliceIndexK);

%Initialize view
pause(1/updateFrequency); %Wait so plot will update
plotSlice([],[],{hf,hSlider_I,1});
hf.UserData.sv3.time=t; %Reset clock so this happens now
plotSlice([],[],{hf,hSlider_J,2}); 
hf.UserData.sv3.time=t; %Reset clock so this happens now
plotSlice([],[],{hf,hSlider_K,3});

drawnow;

%%
varargout{1}=hf;

end

function plotSlice(~,~,inputCell)

hf=inputCell{1};

jSlider=inputCell{2};
dirOpt=inputCell{3};
sliceIndex=round(get(jSlider,'Value'));
hf.UserData.sv3.sliceIndices(dirOpt)=sliceIndex;
sliceIndices=hf.UserData.sv3.sliceIndices;

dt=1/hf.UserData.sv3.updateFrequency; 
t=hf.UserData.sv3.time;
t2=clock;
dtt=etime(t2,t); %Elapsed time

if dtt>dt %If ready to update
    hf.UserData.sv3.time=t2;
    
    M=hf.UserData.sv3.M_plot;
    v=hf.UserData.sv3.v;
    patchType=hf.UserData.sv3.patchTypes{dirOpt};
    
    logicPatch=false(size(M));
    switch dirOpt
        case 1
            logicPatch(sliceIndex,:,:)=1;
        case 2
            logicPatch(:,sliceIndex,:)=1;
        case 3
            logicPatch(:,:,sliceIndex)=1;
    end
    
    figure(hf); %TEMP FIX for bug in MATLAB 2018
    
    if isnan(hf.UserData.sv3.hp(dirOpt))
        [F,V,C]=im2patch(M,logicPatch,patchType);
        [V(:,1),V(:,2),V(:,3)]=im2cart(V(:,2),V(:,1),V(:,3),v);
        V=V+hf.UserData.sv3.origin(ones(size(V,1),1),:);
        hf.UserData.sv3.hp(dirOpt)= gpatch(F,V,C,'none',hf.UserData.sv3.alphaLevel);
    else
        V=get(hf.UserData.sv3.hp(dirOpt),'Vertices');
        switch dirOpt
            case 1
                V(:,2)=(sliceIndex-0.5).*v(1);
            case 2
                V(:,1)=(sliceIndex-0.5).*v(2);
            case 3
                V(:,3)=(sliceIndex-0.5).*v(3);
        end
        set(hf.UserData.sv3.hp(dirOpt),'CData',M(logicPatch)); %Set color data
        set(hf.UserData.sv3.hp(dirOpt),'Vertices',V); %Set vertices
    end
    
    navString=['I: ',num2str(sliceIndices(1)),', J:  ',num2str(sliceIndices(2)),', K: ',num2str(sliceIndices(3))];
    
    set(hf.UserData.sv3.ht,'string',navString);
    
    hf.Name=[hf.UserData.sv3.Name,' ',navString];
end
end

function setScrollSizeFunc(~,~,inputCell)
hf=inputCell{1};
w=inputCell{2};

for q=3:numel(inputCell)
    hSlider=inputCell{q};
    posData=[w*(q-3),0,w,round(hf.Position(4))];
    set(hSlider,'Position',posData);    
end

end

 
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
