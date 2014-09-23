function FEBioPath=getFEBioPath

filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');
try 
    [T]=txtfile2cell(fileName);
    FEBioPath=T{1};
catch
    FEBioPath=[];
end
