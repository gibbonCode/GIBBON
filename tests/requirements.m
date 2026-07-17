function products = requirements(what)
% Wrapper for matlab.codetools.requiredFilesAndProducts,
% to get the required MATLAB products for the GIBBON toolbox.

arguments
    what char {mustBeMember(what,{'tests','demos','help','all'})} = 'all'
end

switch lower(what)
    case 'tests'
        pat = {'tests','TEST*.m'};
    case 'demos'
        pat = {'docs','DEMO*.m'};
    case 'help'
        pat = {'docs','HELP*.m'};
    case 'all'
        p = cellfun(@requirements, {'tests','demos','help'}, 'UniformOutput', false);
        products = unique([p{:}]);
        return;
end

files = cellstr(ls(fullfile(gibbonSettings.gibbonPath, pat{:})));
[~, productStruct] = matlab.codetools.requiredFilesAndProducts(files);
products = {productStruct.Name};
