function [varargout]=plotContours(varargin)
% function [handleCell]=plotContours(contourSet,optionStruct)


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
optionStructDef.hAxis=gca;
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
        loadedData=load(loadName);
        contourData{q}=loadedData.saveStruct.ContourSet;
    end
end

handleCell=cell(numContours,1);
axes(optionStruct.hAxis);
for q=1:1:numContours
    HP=gobjects(1,1);
    c=1;
    Vcs=contourData{q};    
    for qSlice=1:1:numel(Vcs)
        numSubContours=numel(Vcs{qSlice});        
        for qSub=1:1:numSubContours
            Vd=Vcs{qSlice}{qSub}; %Current contour            
            if ~isempty(Vd)                
                hp=plotV(Vd,'w-');
                HP(c,1)=hp;
                c=c+1; 
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
    handleCell{q}=HP;
end

varargout{1}=handleCell;

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
