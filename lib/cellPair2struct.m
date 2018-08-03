function [S]=cellPair2struct(nameCell,entryCell,convertOption)

% Convert cell to structure
S=[];
for q=1:numel(nameCell)
    if convertOption==1
        numberFix = sscanf(entryCell{q},'%f');
    else
        numberFix=[];
    end
    if ~isempty(numberFix)
        S.(nameCell{q})=numberFix;
    else
        S.(nameCell{q})=entryCell{q};
    end
end