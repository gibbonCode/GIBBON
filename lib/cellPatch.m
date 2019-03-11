function HP=cellPatch(varargin)

% % Vv(isinf(Vv))=NaN;
% 
% % L_plot=cellfun(@(X) ~any(isnan(Vv(X(:)))),Fv_cell);
% 
% plotRange=1:numel(Fv_cell);%find(L_plot)
% HP=nan(numel(plotRange),1);
% qI=1;
% for q=plotRange'
%     fv=Fv_cell{q};
%     
%     hp=patch('faces',fv,'vertices',Vv,'faceColor','r');
%     
%         HP(qI)=hp; %Store handle
%     
%         if nargin>2 %Color specified
%             if ischar(Cv) %String input for color
%                 set(hp,'faceColor',Cv);
%             elseif size(Cv,2)==1 %Colormap driven
%                 set(hp,'faceColor','flat','CData',Cv(q,:));
%             elseif size(Cv,2)==3 %RGB driven
%                 set(hp,'faceColor',Cv(q,:));
%             end
%         end
%     qI=qI+1;
% end
% 
% end

switch nargin    
    case 2
        Vv=varargin{1};
        Fv_cell=varargin{2};
        Cv='r';
    case 3
        Vv=varargin{1};
        Fv_cell=varargin{2};
        Cv=varargin{3};
end

Vv(isinf(Vv))=NaN;

L_plot=cellfun(@(X) ~any(isnan(Vv(X(:)))),Fv_cell);

plotRange=find(L_plot);
HP=nan(numel(plotRange),1);
qI=1;
HP=[];
for q=plotRange(:)'
    Ft=Fv_cell{q};
    hp=patch('Faces',Ft,'Vertices',Vv,'FaceColor','b');
    
    if ischar(Cv) %String input for color
        set(hp,'faceColor',Cv);
    elseif size(Cv,2)==1 %Colormap driven
        set(hp,'faceColor','flat','CData',Cv(q,:));
    elseif size(Cv,2)==3 %RGB driven
        set(hp,'faceColor',Cv(q,:));
    end

    HP(qI)=hp;
    qI=qI+1;
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
% Copyright (C) 2019  Kevin Mattheus Moerman
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
