function [varargout]=importFEBio_logfile(varargin)

% function [TIME, DATA, Data_label]=importFEBio_logfile(import_name)
% ------------------------------------------------------------------------
%
% 
% Data format:
%
%   *Title:
%   *Step  = 1
%   *Time  = 0.05
%   *Data  = x;y;z
%   1, -5.01378, 0, 0
%   2, -5.01378, 0, 1.24313
%   3, -5.01378, 0, 2.48626
%   *Step  = 2
%   *Time  = 0.1
%   *Data  = x;y;z
%   1, -5.01378, 0, 0
%   2, -5.01378, 0, 1.24313
%   3, -5.01378, 0, 2.48626
%
% ------------------------------------------------------------------------

%% Parse input

switch nargin
    case 1 
        fileNameImport=varargin{1};
        addZeroInitialOpt=0;
        removeIndicesOpt=0;
    case 2
        fileNameImport=varargin{1};
        addZeroInitialOpt=varargin{2};
        removeIndicesOpt=0;
    case 3
        fileNameImport=varargin{1};
        addZeroInitialOpt=varargin{2};
        removeIndicesOpt=varargin{3};
end

%% Loading .txt file into cell array
[T]=txtfile2cell(fileNameImport);

%% Access data and label parts
logicData=~gcontains(T,'*'); 
T_data=T(logicData);
T_labels=T(~logicData);

%% Getting time increments
logicTimeFlag=gcontains(T_labels,'*Time'); 
% numSteps=nnz(logicTimeFlag);
timeVec=cell2mat(cellfun(@(x) (sscanf(x,'*Time = %f')'),T_labels(logicTimeFlag),'UniformOutput',0));

%% Getting data

logicDataLabel=gcontains(T_labels,'*Data'); 
dataLabel=T_labels{find(logicDataLabel,1)};
dataLabel=regexprep(dataLabel,'*Data  = ','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMP FIX to cope with FEBio bug for Windows users https://github.com/febiosoftware/FEBio/issues/6
if ~isunix
    %Replace '-1.#IND,-1.#IND,-1.#IND' by '-nan,-nan,-nan'
    T_data=regexprep(T_data,'-1.#IND,-1.#IND,-1.#IND','-nan,-nan,-nan');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataArray=cell2mat(cellfun(@(x) (sscanf(x,'%f,')'),T_data,'UniformOutput',0));

if removeIndicesOpt==1
    dataArray=dataArray(:,2:end);
end

if nnz(logicData)>1        
    
    %Check if last export was interupted
    [setSizesAll]=adjacentdircount(logicData,1); %data part sizes   
    setSizes=unique(setSizesAll(setSizesAll>0));     
    numDataEntriesStep=max(setSizes); %The largest (assumed complete) set    
    logicIncomplete=numel(setSizes)>1;
    
    %Trim data to completed parts only if needed
    if logicIncomplete %If more than one size is discovered analysis was likely incomplete        
        dataArray=dataArray(1:end-rem(size(dataArray,1),numDataEntriesStep),:);        
        warning(['Incomplete data detected in: ',fileNameImport,', analysis may have terminated prematurely'])
    end
    
    %Reshape data array into 3D matrix
    dataMat=permute(reshape(permute(dataArray,[2,3,1]),size(dataArray,2),numDataEntriesStep,size(dataArray,1)./numDataEntriesStep),[2,1,3]);
    
    %Trim time vector to match number of entries kept for dataMat
    timeVec=timeVec(1:size(dataMat,3));     
end

%% Add zero initial state if requested

if addZeroInitialOpt==1    
    %Add first time point
    timeVec=[0; timeVec(:)]; 
    
    %Add initial zero data    
    sizImport=size(dataMat);
    sizImport(3)=sizImport(3)+1;
    dataMat_n=zeros(sizImport);
    dataMat_n(:,:,2:end)=dataMat;
    dataMat=dataMat_n;
end

%% Collect output

switch nargout
    case 1
        dataStruct.time=timeVec; 
        dataStruct.data=dataMat;         
        varargout{1}=dataStruct;
    case 2
        varargout{1}=timeVec;
        varargout{2}=dataMat;
        varargout{3}=dataLabel;
    case 3
        varargout{1}=timeVec;
        varargout{2}=dataMat;
        varargout{3}=dataLabel;
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
