function [F,V]=import_obj_geom(fileName)

% function [F,V]=import_obj_geom(fileName)
% -----------------------------------------------------------------------
% This function imports only the geometry data from the obj file defined by
% fileName. 
%
% Change log: 
% 2021/04/30 KMM: Fixed bug relating to / symbols for faces 
% 2021/04/30 KMM: Speeded up by using single cellfun based loop
% -----------------------------------------------------------------------

%% Import file to a cell array
T=txtfile2cell(fileName);

%% Parse cell array for face entries to keep only vertex index 
% i.e. each of these: 
% f v1/vt1 v2/vt2 v3/vt3 ...
% f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 ...
% f v1//vn1 v2//vn2 v3//vn3 ...
% Becomes: 
% v1 v2 v3 ...

T=regexprep(T,'/+[0-9]+|f.', ''); %Remove / signs and anything after up to next space

%% Get faces and vertices

[FC,VC]=cellfun(@(t) parseLineFun(t),T,'UniformOutput',0);
F=cell2mat(FC); %Faces
V=cell2mat(VC); %Vertices

end

%%

function [f,v]=parseLineFun(t)

if numel(t)>2
    switch t(1:2)
        case 'v '
            f=[];
            v=sscanf(t,'v %f %f %f')';
        otherwise
            f=sscanf(t,'%d')';
            v=[];
    end
else
    f=[];
    v=[];
end

end

%% Slower alternative
% 
% T=txtfile2cell(fileName);
% 
% n=numel(T);
% V=nan(n,3);
% c=true(1,1);
% for q=1:1:n
% %     t=strtrim(T{q}); %Text without preceeding/trailing spaces    
% t=T{q};
%     if strcmp(t(1:2),'v ')
%         V(q,:)=sscanf(t,'v %f %f %f');
%     elseif strcmp(t(1:2),'f ')
%         t=regexprep(t,'/+[0-9]+|f.', '');
%         f=sscanf(t,'%d');
%         if c==1
%             F=nan(n,numel(f));
%         end
%         F(q,:)=f;
%     end    
% end
% F=F(~isnan(F(:,1)),:);
% V=V(~isnan(V(:,1)),:);
%%

%% 
% _*GIBBON footer text*_ 
% 
% License: <https://github.com/gibbonCode/GIBBON/blob/master/LICENSE>
% 
% GIBBON: The Geometry and Image-based Bioengineering add-On. A toolbox for
% image segmentation, image-based modeling, meshing, and finite element
% analysis.
% 
% Copyright (C) 2006-2021 Kevin Mattheus Moerman and the GIBBON contributors
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
