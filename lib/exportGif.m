function exportGif(optStruct)

%% Parse input

%Check for files names
if ~isfield(optStruct,'FileNames')    
    error('Please provide FileNames in the input structure');    
end

%Force optStruct.FileNames to be cell
if ~isa(optStruct.FileNames,'cell')
    optStruct.FileNames={optStruct.FileNames};
end

%Check if input is a single folder or list of files
if numel(optStruct.FileNames)==1    
    checkName=optStruct.FileNames{1};    
    if exist(checkName,'dir')==7 
        %It is a folder so attempt to convert all files in folder to a gif.
        %File names are sorted. 
        FileNames=dir(fullfile(checkName));
        FileNames={FileNames(1:end).name};
        FileNames=sort(FileNames(:));
        c=1;
        for q=1:1:numel(FileNames)
            fullfile(checkName,FileNames{q})
            if exist(fullfile(checkName,FileNames{q}),'file')==2
                optStruct.FileNames{c}=fullfile(checkName,FileNames{q});
                c=c+1;
            end
        end
    end
end

%Use defaults if other inputs are missing    
if ~isfield(optStruct,'DelayTime')
    optStruct.DelayTime=0.5; 
end

if ~isfield(optStruct,'LoopCount')
    optStruct.LoopCount=inf; 
end

if ~isfield(optStruct,'FileNameGif')
    optStruct.FileNameGif='exportGif'; 
end

%Check if save location exists
[Savepath,~,fileExt] = fileparts(optStruct.FileNameGif);
if exist(Savepath,'dir')~=7 %create output folder if it does not exist already
    mkdir(Savepath);
end

%Check if extension is given
if isempty(fileExt)    
    optStruct.FileNameGif=[optStruct.FileNameGif,'.gif'];
% elseif ~strcmp(fileExt,'.gif')    
end

%%

numFiles=numel(optStruct.FileNames);
hw = waitbar(0,'Exporting .gif animation');
for q=1:1:numFiles   
    D = imread(optStruct.FileNames{q});    
    [A,map] = rgb2ind(D,256);    
    if q == 1
        imwrite(A,map,optStruct.FileNameGif,'gif','LoopCount',optStruct.LoopCount,'DelayTime',optStruct.DelayTime);
    else
        imwrite(A,map,optStruct.FileNameGif,'gif','WriteMode','append','DelayTime',optStruct.DelayTime);
    end    
    waitbar(q/numFiles,hw,['Exporting .gif animation ',num2str(round(100*q/numFiles)),'%']);
end
close(hw);

 
%% <-- GIBBON footer text --> 
% Copyright 2017 Kevin Mattheus Moerman
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%   http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
