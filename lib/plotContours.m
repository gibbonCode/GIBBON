function plotContours(varargin)
% function [varargout]=plotContours(contourSet,optionStruct)


%% Parse input

switch nargin
    case 1
        contourSet=varargin{1};
        optionStruct=[];        
    case 2
        contourSet=varargin{1};
        optionStruct=varargin{2};
end
if isa(contourSet,'char')
   contourSet={contourSet}; 
end

numContours=numel(contourSet);

optionStructDef.LineWidth=5; 
optionStructDef.Color=gjet(numContours);    
optionStructDef.pathName=[];
[optionStruct]=structComplete(optionStruct,optionStructDef,1); %Complement provided with default if missing

%%

    
if ~isa(contourSet{1},'char')
   contourData=contourSet;
else    
    contourData=cell(numel(contourSet),1);
    for q=1:1:numContours     
        if ~isempty(optionStruct.pathName)
            loadName=fullfile(optionStruct.pathName,contourSet{q});            
        else
            loadName=contourSet{q};
        end        
        load(loadName);        
        contourData{q}=saveStruct.ContourSet;
    end
end

for q=1:1:numContours
    Vcs=contourData{q};
    for qSlice=1:1:numel(Vcs)
        numSubContours=numel(Vcs{qSlice});
        for qSub=1:1:numSubContours
            Vd=Vcs{qSlice}{qSub}; %Current contour            
            if ~isempty(Vd)
                hp=plotV(Vd,'w-');
                if isnumeric(optionStruct.Color)
                    if size(optionStruct.Color,1)>1
                        hp.Color=optionStruct.Color(q,:);
                    else
                        hp.Color=optionStruct.Color;
                    end
                else
                    if numel(optionStruct.Color)>1
                        hp.Color=optionStruct.Color{q};
                    else
                        hp.Color=optionStruct.Color;
                    end
                end
                hp.LineWidth=optionStruct.LineWidth;                
            end
        end
    end
end
