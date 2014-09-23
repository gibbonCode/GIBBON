function pathNameFix=pathSepFix(pathName)

fileSepString=filesep;

switch fileSepString
    case '/' %i.e. PC so fix MAC type
        pathNameFix=regexprep(pathName,'\','/');
    case'\' %i.e. MAC so fix PC type
        pathNameFix=regexprep(pathName,'/','\');        
end