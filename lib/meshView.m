function [varargout]=meshView(varargin)

% function [varargout]=meshView(varargin)
% 
% 2018/01/23 Updated to allow for subfigure plotting

%% Parse input

switch nargin
    case 1        
        meshStruct=varargin{1};
        optionStruct=[];
    case 2
        meshStruct=varargin{1};
        optionStruct=varargin{2};
    otherwise
        error('Wrong number of input arguments');
end

defaultOptionStruct.hFig=[];
defaultOptionStruct.numSLiceSteps=25; %Number of animation steps
defaultOptionStruct.cMap=[];
defaultOptionStruct.faceAlpha1=0.2;
defaultOptionStruct.faceAlpha2=1;
defaultOptionStruct.lightWeightPlot=1;
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,0);

%% Get control parameters

figHandles=optionStruct.hFig;
hFig=figHandles(1);
if numel(figHandles)==2    
    hSub=figHandles(2);
else
    hSub=[];
end

numSliceSteps=optionStruct.numSLiceSteps;
cMap=optionStruct.cMap;
faceAlpha1=optionStruct.faceAlpha1;
faceAlpha2=optionStruct.faceAlpha2;
lightWeightPlot=optionStruct.lightWeightPlot;
fontSize=15;

%% Access mesh data

%                 nodes: [415×3 double]
%         facesBoundary: [320×3 double]
%        boundaryMarker: [320×1 double]
%                 faces: [8144×3 double]
%              elements: [2036×4 double]
%     elementMaterialID: [2036×1 double]
%        faceMaterialID: [8144×1 double]
%        loadNameStruct: [1×1 struct]
   
Fb=meshStruct.facesBoundary;
Cb=meshStruct.boundaryMarker;
V=meshStruct.nodes;
CE=meshStruct.elementMaterialID;
E=meshStruct.elements;
F=meshStruct.faces;
CF=meshStruct.faceMaterialID;

numMaterials=numel(unique(CE(:)));
if isempty(cMap)
   cMap=gjet(numMaterials); 
end
%%
% prepare figure
if isempty(hFig)
    hFig=cFigure; %Store figure handle
    title('Cut view of mesh','FontSize',fontSize);
end

figure(hFig);
if ~isempty(hSub)
    subplot(hSub); 
end
hold on;

if ~isempty(Fb)
    gpatch(Fb,V,0.5*ones(1,3),'none',faceAlpha1);
end

hp=gpatch(Fb,V,Cb,'k',faceAlpha2); %Graphics object to vary property of during animation

camlight headlight;
axisGeom(gca,fontSize); 
if numMaterials>1
    caxis([min(CE(:)) max(CE(:))]);
else
%     caxis([0 1]);
end
colormap(gca,cMap); 
if numMaterials>1
    colorbar;
end

%%
% Set up animation

Y=V(:,2); YE=mean(Y(E),2);

animStruct.Time=linspace(0,1,numSliceSteps); %Time vector
cutLevel=linspace(min(Y(:)),max(Y(:)),numSliceSteps); %Property to set

for q=1:1:numSliceSteps %Step through time       
    cutLevelNow=cutLevel(q); %The current cut level    
    
    logicCutView=YE>cutLevelNow;
    [Fs,Cs]=element2patch(E(logicCutView,:),CE(logicCutView));
    
    if lightWeightPlot==1
        [indBoundary]=tesBoundary(Fs,V);
        Fs=Fs(indBoundary,:);
        Cs=Cs(indBoundary,:);
    end
    
    %Set entries in animation structure
    animStruct.Handles{q}=[hp hp]; %Handles of objects to animate
    animStruct.Props{q}={'Faces','CData'}; %Properties of objects to animate
    animStruct.Set{q}={Fs,Cs}; %Property values for to set in order to animate
end

%Add animation layer
anim8(hFig,animStruct);
set(hFig.UserData.anim8.sliderHandles{1},'Value',round(numSliceSteps/2)); %Set to middle
drawnow;

%%

varargout{1}=hFig;
varargout{2}=hp;

