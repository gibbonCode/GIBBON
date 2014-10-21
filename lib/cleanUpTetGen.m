function cleanUpTetGen(pathNameTempFiles)

extCell={'ele','node','face','edge','mtr','smesh','p2t'}; %Extensions of files to delete

for qc=1:1:numel(extCell)
    ext=extCell{qc}; %Current extension
    fileList = dir(fullfile(pathNameTempFiles,['*.',ext]));
    fileList={fileList(1:end).name}; %Current file list
    
    if ~isempty(fileList)
        %Copying files to output location
        for q=1:1:numel(fileList);
            fileName=fullfile(pathNameTempFiles,fileList{q});
            delete(fileName);
        end
    end
end
