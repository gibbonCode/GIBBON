function setFEBioPath(FEBioPathSpec)

%%
filePath=mfilename('fullpath');
toolboxPath=fileparts(fileparts(filePath));
configPath=fullfile(toolboxPath,'config');

fileName=fullfile(configPath,'FEBioPath.txt');

%Import text file containing paths

T{1}='#List paths for febio here, multiple can be given e.g. for use on several platforms. The first valid path found is always used.';

%Set paths in cell
switch class(FEBioPathSpec)
    case 'cell'
        for q=1:1:numel(FEBioPathSpec)
            pathNameNow=FEBioPathSpec{q};
            if exist(pathNameNow,'file')==2
                T{q+1}=FEBioPathSpec{q};
            end
        end
    otherwise
        if exist(FEBioPathSpec,'file')==2
            T{2}=FEBioPathSpec;
        end
end

cell2txtfile(fileName,T,0); %Write cell to config file

end


