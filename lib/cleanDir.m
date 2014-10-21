function cleanDir(pathName)

fileList = dir(pathName); 
fileList={fileList(1:end).name}; 

for q=1:1:numel(fileList);
    fName=fullfile(pathName,fileList{q}); 
    if ~isdir(fName)        
        delete(fName);        
    end
end
