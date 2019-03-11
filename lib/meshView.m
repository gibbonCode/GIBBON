function [varargout]=meshView(varargin)

% function [varargout]=meshView(varargin)
% ------------------------------------------------------------------------
% 2018/01/23 Updated to allow for subfigure plotting
% 2018/04/16 Added control of direction of cutting (X,Y, or Z) and also
% what side is viewed/cut away. 
% 2019/09/07 Fixed bug in relation to single element/face input
% 2019/09/07 Updated handling of cutting direction
% ------------------------------------------------------------------------
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
defaultOptionStruct.numSliceSteps=25; %Number of animation steps
defaultOptionStruct.cMap=[];
defaultOptionStruct.cutDir=2;
defaultOptionStruct.cutSide=1;
defaultOptionStruct.faceAlpha1=0.2;
defaultOptionStruct.faceAlpha2=1;
defaultOptionStruct.lightWeightPlot=1;
[optionStruct]=structComplete(optionStruct,defaultOptionStruct,0);

%% Get control parameters

figHandles=optionStruct.hFig;
if ~isempty(figHandles)
    hFig=figHandles(1);
    if numel(figHandles)==2
        hSub=figHandles(2);
    else
        hSub=[];
    end
else
    hFig=[];
    hSub=[];
end

numSliceSteps=optionStruct.numSliceSteps;
cMap=optionStruct.cMap;
cutDir=optionStruct.cutDir; 
cutSide=optionStruct.cutSide; 
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

E=meshStruct.elements;
V=meshStruct.nodes;

if isfield(meshStruct,'facesBoundary')
    Fb=meshStruct.facesBoundary;
    Cb=meshStruct.boundaryMarker;
else
    F=meshStruct.faces;
    indBoundary=tesBoundary(F,V);
    Fb=F(indBoundary,:);
    Cb=ones(size(Fb,1),1);
end

if isfield(meshStruct,'elementMaterialID')
    CE=meshStruct.elementMaterialID;
else
    CE=ones(size(E,1),1);
end

% F=meshStruct.faces;
% CF=meshStruct.faceMaterialID;

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
    icolorbar;
end

%%
% Set up animation

[VE]=patchCentre(E,V);
XE=VE(:,cutDir);
animStruct.Time=linspace(0,1,numSliceSteps); %Time vector
cutLevel=linspace(min(XE)-max(eps(XE)),max(XE)+max(eps(XE)),numSliceSteps); %Property to set

for q=1:1:numSliceSteps %Step through time       
    cutLevelNow=cutLevel(q); %The current cut level    
    
    if cutSide==1
        logicCutView=XE>=cutLevelNow;
    elseif cutSide==-1
        logicCutView=XE<=cutLevelNow;
    end
    
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

%% Collect output

switch nargout
    case 1
        varargout{1}=hFig;
    case 2
        varargout{1}=hFig;
        varargout{2}=hp;
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
