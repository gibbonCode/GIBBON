function tripletFileRename(folderName,fileExtensions,tripletNamekeepIndices,tripletNameEndIndex,tripletNameCombineIndex,runOption)

%%REMOVE COMMENTING FOR DEBUGGING:
% clear all; close all; clc;
%
% fileExtensions={'.PAR','.REC'};
% folderName='C:\Users\kmmoerman\Desktop\temp';
% tripletNamekeepIndices=[1:5 7]; %triplet combined name selector
% tripletNameEndIndex=6; %This is where p, m or s is found
% tripletNameCombineIndex=8; %This is which column will form the combined filename

for qq=1:1:numel(fileExtensions)
    fileExtension=fileExtensions{qq};
    fileNameCell = dir(fullfile(folderName,['*',fileExtension]));
    fileNameCell={fileNameCell(1:end).name};
    splitFileNames=regexp(fileNameCell, '_', 'split')';
    splitFileNamesMat=reshape([splitFileNames{:}],numel(splitFileNames{1}),numel(splitFileNames))';
    splitFileNamesCellDoubles=cellfun(@double,splitFileNamesMat,'UniformOutput',0);
    splitFileNamesMatDoubles=str2double(splitFileNamesMat);
    
    %% Sorting parts to combine
    
    toSort=splitFileNamesMatDoubles(:,tripletNameCombineIndex);
    [~,sortInd]=sort(toSort);
    fileNameCell=fileNameCell(sortInd);
    splitFileNames=splitFileNames(sortInd);
    splitFileNamesMat=splitFileNamesMat(sortInd,:);
    
    
    %% Finding common name elements
    
    % %Convert to cell array containing doubles
    % stringDoubleReps=cell(size(splitFileNames));
    % for q1=1:1:numel(splitFileNames)
    %    for q2=1:1:numel(splitFileNames{q1})
    %         stringDoubleReps{q1}{q2}=double(splitFileNames{q1}{q2});
    %    end
    % end
    %
    % %N.B. assuming the same number of split components for triplet entries
    % q1=1;%assuming that if difference exists that it exists for all entries
    % L=zeros(1,numel(splitFileNames{q1}));
    % for q2=1:1:numel(splitFileNames{q1})
    %     a=stringDoubleReps{q1}{q2}; b=stringDoubleReps{q1+1}{q2};
    %     L(q2)=any(a~=b); %Are any parts of the split string entries different
    % end
    % differenceIndex=find(L);
    
    %%
    
    %Creating string for beginning of file name, assuming choice is same for
    %each entry in triplet and only first is used to generate the string here
    tripletNamePartCustum=[];
    for q2=tripletNamekeepIndices
        tripletNamePartCustum=[tripletNamePartCustum,splitFileNames{1}{q2},'_'];
    end
    tripletNamePartCustum=tripletNamePartCustum(1:end-1);
    
    %Creating combined file name part
    tripletNamePartCombine=[];
    for q3=1:3
        tripletNamePartCombine=[tripletNamePartCombine,splitFileNamesMat{q3,tripletNameCombineIndex},'_'];
    end
    tripletNamePartCombine=tripletNamePartCombine(1:end-1);
    
    %Creating total new file name
    fileNameCellNew=fileNameCell;
    for q1=1:1:numel(splitFileNames)
        fileNameCellNew{q1}=[tripletNamePartCustum,'_',tripletNamePartCombine,'_',splitFileNames{q1}{tripletNameEndIndex},fileExtension];
    end
    
    %% Renaming the files
    for q2=1:1:numel(fileNameCell);
        oldFileName=fullfile(folderName,fileNameCell{q2}); %The old file name
        newFileName=fullfile(folderName,fileNameCellNew{q2}); %The new file name
        if runOption==1            
                movefile(oldFileName,newFileName); %Overwriting files with their new name           
        end
        disp([fileNameCell{q2},' ----> ',fileNameCellNew{q2}]);
    end
    
end