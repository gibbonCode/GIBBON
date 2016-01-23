function changeFileExtensions(pathName,extOld,extNew)

%%
if isempty(extOld)
    files = dir(pathName);
    files={files(1:end).name};
else
    files = dir(fullfile(pathName,['*.',extOld]));
    files={files(1:end).name};
end

for q=1:1:numel(files)
    oldNameFull=fullfile(pathName,files{q});
    if ~isdir(oldNameFull)
        [~,oldName,extFile] = fileparts(oldNameFull);
        if isempty(extNew)
            newNameFull=fullfile(pathName,oldName);
        else
            newNameFull=fullfile(pathName,[oldName,'.',extNew]);
        end
        if strcmp(oldNameFull,newNameFull)==0
            if isempty(extOld)
                if strcmp(extFile,extOld);
                    movefileNow(oldNameFull,newNameFull);
                end
            else
                movefileNow(oldNameFull,newNameFull);
            end
        end
    end
end

end

function movefileNow(oldNameFull,newNameFull)
    movefile(oldNameFull,newNameFull);
end