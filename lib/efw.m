function efw(varargin)

% function efw(hf)
% ------------------------------------------------------------------------
% Export Figure Widget
%
% Kevin Mattheus Moerman
% gibbon.toolbox@gmail.com
% 
%------------------------------------------------------------------------

%% Parse input arguments
switch nargin
    case 0
        hf = gcf;
    case 1
        hf=varargin{1};
    otherwise
        error('Wrong number of input arguments');
end

if ~ishandle(hf)
    hf = gcf;
end


%% Initialise button
hb = findall(hf,'Type','uitoolbar');

%Check for presence of a efw button
hp = findobj(hb,'Tag','efw_button');
if isempty(hp); %If efw button is not present create one 
    
    % Build icon
    s=[NaN,NaN,0.02,0.64,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.64,0.01,NaN,NaN;...
       NaN,NaN,0.4,0.7,0.15,0.16,0.15,0.15,0.16,0.16,0.16,0.15,0.71,0.38,NaN,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.58,0.43,NaN,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.58,0.42,NaN,NaN;...
       NaN,NaN,0.45,0.58,0.03,0.05,0.05,0.05,0.05,0.05,0.05,0.03,0.6,0.44,NaN,NaN;...
       0.33,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.33;...
       0.78,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.48,0.16,1.0,0.78;...
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.38,0.01,0.79,1.0;...
       1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0;...
       1.0,1.0,0.78,0.62,0.16,0.17,0.17,0.17,0.17,0.17,0.17,0.16,0.63,0.79,1.0,1.0;...
       0.79,1.0,0.78,0.56,NaN,NaN,1,1,1,1,NaN,NaN,0.58,0.78,1.0,0.79;...
       0.51,1.0,0.78,0.56,NaN,NaN,1,NaN,NaN,1,NaN,NaN,0.58,0.79,1.0,0.51;...
       NaN,0.18,0.54,0.56,NaN,NaN,1,1,1,1,NaN,NaN,0.58,0.53,0.18,NaN;...
       NaN,NaN,0.45,0.56,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN,0.58,0.43,NaN,NaN;...
       NaN,NaN,0.4,0.7,0.15,0.15,1,1,1,1,0.16,0.15,0.71,0.38,NaN,NaN;...
       NaN,NaN,0.02,0.64,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.63,0.02,NaN,NaN];
    
    S=zeros(16,16,3);
    S(:,:,1)=0.7.*s;
    S(:,:,2)=0.25.*s;
    S(:,:,3)=0.05.*s;
    
    % Create a uipushtool in the toolbar
    hp=uipushtool(hb,'TooltipString','Export Figure Widget','CData',S,'Tag','efw_button');
    
    defStruct=get(hb,'UserData');
    if isempty(defStruct)
        defStruct.defaultPath=cd;
        defStruct.imName=['fig',num2str(get(hf,'Number'))];
        defStruct.imExt='png';
        defStruct.imRes='120';
        defStruct.exportFigOpt='';
        set(hb,'UserData',defStruct);
    end
    set(hp,'ClickedCallback',{@start_efw_push,{hf,hb}});  
    
end

return

function start_efw_push(hObject,callbackdata,inputData)
start_efw(inputData{1},inputData{2});
return

function start_efw(hf,hb)
figure(hf);
defStruct=get(hb,'UserData');
prompt = {'Save path (leave empty to browse to desired folder instead):','Image name:','Image extension (i.e. png, jpg, pdf, eps, bmp or tif):','Image resolution (e.g. 120):','Extra export_fig options (comma seperated, no spaces e.g. -nocrop,-transparent):'};
dlg_title = 'Export Figure Widget (see: help efw and help export_fig)';
defaultOptions = {defStruct.defaultPath,defStruct.imName,defStruct.imExt,defStruct.imRes,defStruct.exportFigOpt};
s=25+max([cellfun(@numel,prompt) cellfun(@numel,defaultOptions)]);

Q = inputdlg(prompt,dlg_title,[1 s],defaultOptions);

if ~isempty(Q)
    if isempty(Q{1})
        Q{1}=uigetdir(defStruct.defaultPath,'Select save path');
        if Q{1}==0
            Q{1}=[];
        end
    end
    if all(~cellfun(@isempty,Q(1:end-1)))
        fileName=fullfile(Q{1},[Q{2},'.',Q{3}]);
        figRes=['-r',Q{4}];
        inputCell{1,1}=fileName;
        inputCell{1,2}=figRes;
        if ~isempty(Q{5})
            inputCellExtra = strsplit(Q{5},',');
            for q=1:1:numel(inputCellExtra)
                inputCell{end+1}=inputCellExtra{q};
            end
        end
        export_fig(inputCell{:});
        
        defStruct.defaultPath=Q{1};
        defStruct.imName=Q{2};
        defStruct.imExt=Q{3};
        defStruct.imRes=Q{4};
        defStruct.exportFigOpt=Q{5};
        set(hb,'UserData',defStruct);
    else
        return
    end
end

return
