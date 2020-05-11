function gdrawnow
% function gdrawnow
% ------------------------------------------------------------------------
% The |gdrawnow| function is similar to the |drawnow| command but also
% activates the |vcw| widget if present for the current figure window. 
% 
%
% Kevin Mattheus Moerman
%
% Change log: 
% 2020/05/11 Created
% 2020/05/11 Added triple call to drawnow which appears to fix MATLAB bug
% when multiple figures are used and plotting occurs in the wrong figure
% window. 
% 2020/05/11 Added new axis colorbar handling
% ------------------------------------------------------------------------
%%
% Force full drawnow
drawnow; drawnow;  drawnow; %Triple to avoid MATLAB bug

%% Activate vcw 

hf=gcf;

if isfield(hf.UserData,'vcw')
    
    switch hf.UserData.vcw.buttonHandle.State
        case 'on'
            H=findobj(hf,'Type','colorbar'); %Handle set for all colorbars in figure;
            if ~isempty(H)                
                if ~isempty(hf.UserData.vcw.colorbarHandles)
                    Hp=hf.UserData.vcw.colorbarHandles;                    
                    logicSetManual=~ismember(H,Hp);
                    if nnz(logicSetManual)>0
                        H_new=H(logicSetManual);
                        hf.UserData.vcw.colorbarHandles=[hf.UserData.vcw.colorbarHandles H_new];
                        hf.UserData.vcw.colorbarLocSet(end:end+numel(H_new))={H_new.Location};
                        for q=1:1:numel(H_new) %Loop over colorbar handles and set location property
                            set(H_new(q),'Location','manual');
                        end
                    end
                end
            end
        case 'off'
            hf.UserData.vcw.buttonHandle.State='on';
    end
   
    %Set default views for new axes and attempt to hide interactive toolbars 
    h = findobj(hf, 'Type', 'axes', '-depth', 1)'; %All axis handles
    if ~isempty(h)
        for hNow = h
            %Set default views if not present                  
            if ~isfield(hNow.UserData,'defaultView')
                hNow.UserData.defaultView=camview(hNow);
            end            
            
            %Try to hide toolbar
            try
                hNow.Toolbar.Visible = 'off';
            catch
            end
        end
    end
    
    
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
% Copyright (C) 2006-2020 Kevin Mattheus Moerman
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
