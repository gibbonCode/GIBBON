function [varargout]=febioOutputView(pathName,febio_spec)

% function [hf]=febioOutputView(pathName,febio_spec)
% ------------------------------------------------------------------------
%
%
%
% ------------------------------------------------------------------------

%%
% optionStruct=[];

%%
fontSize=15;

%%
hf=cFigure; %Open figure

%%

V=febio_spec.Geometry.Nodes{1}.node.VAL;
E_cell=febio_spec.Geometry.Elements;

%%

nodeDataCellAll=febio_spec.Output.logfile.node_data;

nodalDisplacementAll=zeros(size(V));
for q=1:1:numel(nodeDataCellAll)
    dataTypeString=nodeDataCellAll{q}.ATTR.data;
    if strcmp(dataTypeString,'ux;uy;uz')
        [~,nodalDisplacementAll,~]=importFEBio_logfile(fullfile(pathName,nodeDataCellAll{q}.ATTR.file)); %Import data
        
        nodalDisplacementAll=nodalDisplacementAll(:,2:end,:);
        sizImport=size(nodalDisplacementAll);
        sizImport(3)=sizImport(3)+1;
        dataArray_n=zeros(sizImport);
        dataArray_n(:,:,2:end)=nodalDisplacementAll;
        nodalDisplacementAll=dataArray_n;
        U_end=nodalDisplacementAll(:,:,end);
        U_mag=sqrt(sum(U_end(:,3).^2,2));
        V_def=V+U_end;       
    end
end
nodalDisplacementMagnitudeAll=sqrt(sum(nodalDisplacementAll.^2,2));

hf.UserData.febioOutputView.colorData=[];
hf.UserData.febioOutputView.pathName=pathName;
hf.UserData.febioOutputView.febio_spec=febio_spec;

%%
% Create basic view and store graphics handle to initiate animation

hf.UserData.febioOutputView.ht=gtitle(' - '); %Figure axis title

handleSet=gobjects(1,numel(E_cell));
for q=1:1:numel(E_cell)
    
    elementsNow=febio_spec.Geometry.Elements{q}.elem.VAL;
    elementTypeNow=febio_spec.Geometry.Elements{1}.ATTR.type;
    
    if any(strcmp(elementTypeNow,{'tri3','tri6','tri7','quad4','quad8'})) %Surface elements
        facesNow=elementsNow;
    else %Solid elements
        facesNow=element2patch(elementsNow);
    end
    
    handleSet(q)=gpatch(facesNow,V,U_mag,'k',1); %Add graphics object to animate
    handleSet(q).FaceColor='interp';
end

axisGeom(gca,fontSize);
colormap(gjet(250)); colorbar;
cLim=[min(nodalDisplacementMagnitudeAll(:)) max(nodalDisplacementMagnitudeAll(:))];
if abs(cLim(2)-cLim(1))<eps(1)
    cLim(1)=cLim(1)-1;
    cLim(2)=cLim(2)+1;
end
hf.UserData.febioOutputView.cLim=cLim;
caxis(cLim);
VV=[V;V_def];
axis([min(VV(:,1)) max(VV(:,1)) min(VV(:,2)) max(VV(:,2)) min(VV(:,3)) max(VV(:,3))]); %Set axis limits statically
camlight headlight;

hf.UserData.febioOutputView.handleSet=handleSet;


%%
        
outputDataCellAll=febio_spec.Output.logfile.node_data;
if isfield(febio_spec.Output.logfile,'element_data')
    outputDataCellAll=[outputDataCellAll febio_spec.Output.logfile.element_data];
end


popUpString=cell(1,1);
%Initialize with displacement magnitude
popUpString{1}='um';
hf.UserData.febioOutputView.data.um=nodalDisplacementMagnitudeAll;

for q=1:1:numel(outputDataCellAll)
    
    dataTypeString=outputDataCellAll{q}.ATTR.data;
    dataTypeSet=strsplit(dataTypeString,';');
    
    [timeVec,dataArray,~]=importFEBio_logfile(fullfile(pathName,outputDataCellAll{q}.ATTR.file)); %Import data
    dataArray=dataArray(:,2:end,:);
    sizImport=size(dataArray);
    sizImport(3)=sizImport(3)+1;
    dataArray_n=zeros(sizImport);
    dataArray_n(:,:,2:end)=dataArray;
    dataArray=dataArray_n;
    
    for qSet=1:1:numel(dataTypeSet)
        dataTypeStringNow=dataTypeSet{qSet};
        popUpString{end+1}=dataTypeStringNow;
        dataNow=dataArray(:,qSet,:);
        switch dataTypeStringNow
            case 'J' %Volume ratio
                dataNow(:,:,1)=1; %set initial to ones
        end
        hf.UserData.febioOutputView.data.(dataTypeStringNow)=dataNow;
    end    

end

timeVec=[0; timeVec(:)]; %Time

%%


backGroundColor=hf.Color;
if strcmp(backGroundColor,'k')
    textColor='w';
elseif size(backGroundColor,2)==3
    grayLevels=linspace(0,1,100);
    [~,indMax]=max(abs(grayLevels-mean(double(backGroundColor))));
    textColor=grayLevels(indMax)*ones(1,3);
else
    textColor='k';
end

hPop = uicontrol(hf,'Style','popupmenu','String',popUpString,...
    'BackgroundColor',backGroundColor,'FontSize',fontSize,'FontWeight','normal','Units','Points','ForegroundColor',textColor,...
    'Callback',@(a,b)setColorData(a,b,{hf}));

hFunc=get(hf,'ResizeFcn');

if iscell(hFunc)
    warning('gtitle replaced the ResizeFcn function. Specify your ResizeFcn in the form @(h,e)figResize(h,e,c) to avoid this behavior');    
    set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hPop}));
else
    if isempty(hFunc)
        set(hf,'ResizeFcn',@(a,b)figResize(a,b,{hf,hPop}));
    else        
        set(hf,'ResizeFcn',@(a,b)(cellfun(@(x)feval(x,a,b),{hFunc,@(a,b)figResize(a,b,{hf,hPop})})));
    end
end

figResize([],[],{hf,hPop});

%%

clear animStruct;

% Set up animation features
animStruct.Time=timeVec; %The time vector
propertySetVert=repmat({'Vertices'},1,numel(handleSet));
propertySetCData=repmat({'CData'},1,numel(handleSet));

for qt=1:1:size(nodalDisplacementAll,3) %Loop over time increments
    U_now=nodalDisplacementAll(:,:,qt); %Current displacement
    V_def=V+U_now; %Current nodal coordinates    
    d1=repmat({V_def},1,numel(handleSet));
    
    %Initialize animation structure with nodal coordinate and displacement animation
    animStruct.Handles{qt}=[handleSet]; %Handles of objects to animate
    animStruct.Props{qt}={propertySetVert{:}}; %Properties of objects to animate
    animStruct.Set{qt}={d1{:}}; %Property values for to set in order to animate
      
    %Add colordata animations for faces
    colorData=nodalDisplacementMagnitudeAll(:,:,qt);
    
    d2={};
    for qs=1:1:numel(handleSet)
        colorDataNow=colorData;        
%         if size(colorData,1)==size(V,1)            
%             f=handleSet(qs).Faces;
%             [colorDataNow]=vertexToFaceMeasure(f,colorData);
%         end
        f=handleSet(qs).Faces;
        if size(colorData,1)==size(f,1)                    
            [colorDataNow]=faceToVertexMeasure(f,V,colorData);            
        end
        d2{qs}=colorDataNow;        
    end    
    animStruct.Handles{qt}(1:2*numel(handleSet))=[handleSet handleSet]; %Handles of objects to animate
    animStruct.Props{qt}(end+1:end+numel(propertySetCData))={propertySetCData{:}}; %Properties of objects to animate
    animStruct.Set{qt}(end+1:end+numel(propertySetCData))={d2{:}}; %Property values for to set in order to animate
end

anim8(hf,animStruct); %Initiate animation feature
drawnow;

%%

if nargout>0
    varargout{1}=hf;
end

end


%%

function figResize(~,~,inputCell)

hf=inputCell{1};
hPop=inputCell{2};
unitsNow=hf.Units;
hf.Units='Points';

figPosition=hf.Position;
textBoxWidth=round(figPosition(3)/6);
textBoxHeight=hPop.Extent(4).*ceil(hPop.Extent(3)./figPosition(3));
textPosition=[0 figPosition(4)-textBoxHeight textBoxWidth textBoxHeight];
hPop.Position = textPosition;

hf.Units=unitsNow;

end

%%

function setColorData(source,~,inputCell)
hf=inputCell{1};
indexValue = source.Value;
stringSet = source.String;

dataTypeStringNow = stringSet{indexValue};
hf.UserData.febioOutputView.ht.String=dataTypeStringNow;

hf.UserData.febioOutputView.colorData=hf.UserData.febioOutputView.data.(dataTypeStringNow);

updateFaceColor(hf);

caxis(hf.UserData.febioOutputView.cLim);
drawnow;

end

%%

function updateFaceColor(hf)

for qt=1:1:numel(hf.UserData.anim8.animStruct.Time)
    cLim=[];    
    for qh=numel(hf.UserData.febioOutputView.handleSet)+1:1:2*numel(hf.UserData.febioOutputView.handleSet)        
        c=hf.UserData.febioOutputView.colorData(:,:,qt);
        d2=repmat({c},1,numel(hf.UserData.febioOutputView.handleSet));
        hf.UserData.anim8.animStruct.Set{qt}{qh}=d2{:}; %Property values for to set in order to animate
        
        if isempty(cLim)
            cLim=[min(c(:)) max(c(:))];
        else
            cLim(1)=min(min(c(:)),cLim(1));
            cLim(2)=max(max(c(:)),cLim(2));
        end
    end
end

if abs(cLim(2)-cLim(1))<eps(1)
    cLim(1)=cLim(1)-1;
    cLim(2)=cLim(2)+1;
end
hf.UserData.febioOutputView.cLim=cLim;

%Briefly toggle slider to trigger anim8 redraw event
jSlider=hf.UserData.anim8.sliderHandles{1};
v=get(jSlider,'Value');
if v>1
    set(jSlider,'Value',v-1);
    set(jSlider,'Value',v);
else
    set(jSlider,'Value',v+1);
    set(jSlider,'Value',v);
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
