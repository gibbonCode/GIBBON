function pathNameFix=pathSepFix(pathName)

fileSepString=filesep;

switch fileSepString
    case '/' %i.e. PC so fix to MAC type
        pathNameFix=regexprep(pathName,'\','/');
    case'\' %i.e. LINUX or MAC so fix to PC type
        pathNameFix=regexprep(pathName,'/','\');        
end