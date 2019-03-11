function exportGifAnim8(varargin)


%%

switch nargin
    case 0
        hf=gcf;
        inputStruct=[];
        dlgOpt=0;
    case 1
        hf=varargin{1};
        inputStruct=[];
        dlgOpt=0;
    case 2
        hf=varargin{1};
        inputStruct=varargin{2};
        dlgOpt=0;
    case 3
        hf=varargin{1};
        inputStruct=varargin{2};
        dlgOpt=varargin{3};
end

%%
        
defStruct=hf.UserData.efw;

[inputStruct]=structComplete(inputStruct,defStruct,1); %Complement provided with default if missing or empty

if dlgOpt==1
    prompt = {'Save path (leave empty to browse to desired folder instead):',...
        'Image name:','Image extension (i.e. png, jpg, bmp, or tif):',...
        'Image resolution (e.g. 120):',...
        'Extra export_fig options (comma seperated, no spaces e.g. -nocrop,-transparent,-painters):',...
        'Export gif option:'};
    dlg_title = 'Export Gif Widget (see: help efw and help export_fig)';
    defaultOptions = {inputStruct.defaultPath,...
        inputStruct.imName,...
        inputStruct.imExt,...
        inputStruct.imRes,...
        inputStruct.exportFigOpt,...
        hf.UserData.efw.exportGifOpt,...
        };
    
    s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);
    
    Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);
else
    Q = {inputStruct.defaultPath,...
        inputStruct.imName,...
        inputStruct.imExt,...
        inputStruct.imRes,...
        inputStruct.exportFigOpt,...
        '1'};
end

if ~isempty(Q)
    if isempty(Q{1})
        Q{1}=uigetdir(inputStruct.defaultPath,'Select save path');
        if Q{1}==0
            return;
        end
    end
    
    if isempty(Q{2})
        error('Empty input. Please enter a file name');
    end
    
    if ~exist(Q{1},'dir') %create output folder if it does not exist already
        mkdir(Q{1});
    end
    
    if isempty(Q{3})
        warning('Empty input. No image format specified, using jpg files.');
        Q{3}='jpg';
    end
    
    fileName=fullfile(Q{1},Q{2});
    exportGifCell{1,1}=fileName;
    
    stringSet=Q{3}; %The image extension
    stringNoSpaces=regexprep(stringSet,'[^\w'']',''); %Remove potential extra spaces
    
    if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'
        stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
    end
    
    %Check format validaty and keep if valid
    if any(strcmp(stringNoSpaces,{'-png','-jpg','-tiff','-bmp'}))
        exportGifCell{1,end+1}=stringNoSpaces; %Add to input list
    else
        error('Wrong image format requested');
    end
    
    figRes=['-r',Q{4}];
    exportGifCell{1,end+1}=figRes;
    
    if ~isempty(Q{5})
        stringSet=Q{5}; %The set of potentially multiple options
        stringSetSep = strsplit(stringSet,',');
        for q=1:1:numel(stringSetSep)
            stringNoSpaces=regexprep(stringSetSep{q},'[^\w'']',''); %Remove potential extra spaces
            if ~strcmp(stringNoSpaces(1),'-') %If first character is not '-'
                stringNoSpaces=['-',stringNoSpaces]; %Add '-' to start, e.g. 'jpg' becomes '-jpg'
            end
            exportGifCell{1,end+1}=stringNoSpaces; %Add to input list
        end
    end
    
    if ~isempty(Q{6})
        exportGifOpt=Q{6};
    end
    
    fileNameGif=exportGifCell{1,1};
    exportGifCellSub=exportGifCell;
    
    c=1;
    stepRange=1:hf.UserData.anim8.shiftMag:numel(hf.UserData.anim8.animStruct.Time);    
    for q=stepRange
        set(hf.UserData.anim8.sliderHandles{1},'Value',q);
        fileNameNow=[fileNameGif,'_',num2str(q)];
        exportGifCellSub{1,1}=fileNameNow;
        figure(hf);
        export_fig(exportGifCellSub{:});
        gifStruct.FileNames{c}=[fileNameNow,'.',exportGifCell{1,2}(2:end)];
        c=c+1;
    end
    
    if strcmp(exportGifOpt,'1')
        %Add reverse path
        if strcmp(get(hf.UserData.anim8.ButtonHandles.hCycle,'State'),'on')
            numFiles=numel(gifStruct.FileNames);
            if numFiles>2
                for q=(numFiles-1):-1:2
                    gifStruct.FileNames{end+1}=gifStruct.FileNames{q};
                end
            end
        end
        
        gifStruct.DelayTime=hf.UserData.anim8.pauseTime;
        gifStruct.FileNameGif=fileNameGif;
        
        exportGif(gifStruct);
        
        %Cleanup image files
        for q=1:1:numel(gifStruct.FileNames)
            if exist(gifStruct.FileNames{q},'file')==2
                delete(gifStruct.FileNames{q});
            end
        end
    end
    
    %Override defaults
    defStruct.defaultPath=Q{1};
    defStruct.imName=Q{2};
    defStruct.imExt=Q{3};
    defStruct.imRes=Q{4};
    defStruct.exportFigOpt=Q{5};
    defStruct.efw.exportGifOpt=Q{6};
    
    hf.UserData.efw=defStruct;
    
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
