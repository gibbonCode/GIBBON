function my_builddocsearchdb(helploc)
%BUILDDOCSEARCHDB Build documentation search database
%    BUILDDOCSEARCHDB HELP_LOCATION builds a search database for MYTOOLBOX 
%    documentation in HELP_LOCATION, which the Help browser uses to perform 
%    searches that include the documentation for My Toolbox. Use 
%    BUILDDOCSEARCHDB for My Toolbox HTML help files you add to the Help 
%    browser via an INFO.XML file. BUILDDOCSEARCH creates a directory named
%    HELPSEARCH in HELP_LOCATION. The HELPSEARCH directory works only with 
%    the version of MATLAB used to create it.
%
%    Examples:
%    builddocsearchdb([matlabroot '/toolbox/mytoolbox/help']) - builds the
%    search database for the documentation found in the directory
%    /toolbox/mytoolbox/help under the MATLAB root.
%
%    builddocsearchdb D:\Work\mytoolbox\help - builds the search database
%    for the documentation files found at D:\Work\mytoolbox\help.

%   Copyright 2006-2012 The MathWorks, Inc.

if ~usejava('jvm')
    error('MATLAB:doc:CannotBuildSearchDb','%s',getString(message('MATLAB:doc:UnsupportedPlatform',upper(mfilename))));
end
if (~exist(helploc,'file'))
    error('MATLAB:doc:CannotBuildSearchDb','%s',getString(message('MATLAB:doc:SpecifiedDirectoryDoesNotExist')));
end
try
    com.mathworks.mlwidgets.help.search.lucene.LuceneIndexServices.indexDoc(helploc);
catch 
    error('MATLAB:doc:CannotBuildSearchDb','%s',getString(message('MATLAB:doc:CouldNotWriteSearchDatabase')));
end
if (~com.mathworks.mlwidgets.help.HelpInfo.isInstalledProductHelpLocation(helploc))
    warning(message('MATLAB:doc:DocNotInstalled'));
end
disp(getString(message('MATLAB:doc:CreatedSearchDatabase')));