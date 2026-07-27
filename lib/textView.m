function [varargout]=textView(T, optionStruct, addOpts)

%%
arguments
    T
    optionStruct (1,1) struct = struct()
    addOpts.BackgroundColor (1,3) double = [44 48 55]/255;
    addOpts.FontColor (1,3) double = [168 177 190]/255;
    addOpts.HorizontalAlignment char = 'Left';
    addOpts.FontName char = 'Monospaced';
    addOpts.FontWeight char = 'Bold';
    addOpts.FontSize (1,1) double = 12;
    addOpts.WordWrap char = 'off';
    addOpts.Editable char = 'off';
end

fileName='';
if ~iscell(T)
    if exist(T,'file')==2
        fileName=T;
        T=txtfile2cell(T);
    end
end

% Accept ForegroundColor and option structure, for backwards compatibility
if isfield(optionStruct,'ForegroundColor')
    optionStruct = renameStructField(optionStruct,'ForegroundColor', 'FontColor');
end
optionStruct = structComplete(optionStruct, addOpts, 1);

%%

% Add line numbers
t=(1:1:numel(T));
n=numel(sprintf('%d',max(t)));
t=sprintf(['%-',num2str(n),'.0d    \n'],t); 
t=t(1:end-1);
t=strsplit(t,'\n')';
T=strcat(t,T);

% Open figure
figStruct.Name=fileName;
figStruct.Color=optionStruct.BackgroundColor;
figStruct.MenuBar='none';
figStruct.vcw=0;
figStruct.efw=0;
hf = cFigure(figStruct);

hGrid = uigridlayout(hf, [1 1]);
hGrid.Padding = [0 0 0 0];
hGrid.RowSpacing = 0;
hGrid.ColumnSpacing = 0;

optArgs = [fieldnames(optionStruct), struct2cell(optionStruct)]';
uitextarea(hGrid, 'Value', T, optArgs{:});
drawnow;

%%
if nargout>0
    varargout{1}=hf;
end

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2026 Kevin Mattheus Moerman and the GIBBON contributors
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
