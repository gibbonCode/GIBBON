function [C,group_sizes]=contour_group(varargin)

%function [C,group_sizes]=contour_group(X,Y,Z,M,Sx,Sy,Sz,c,interpMethod)

%% Parse input

switch nargin
    case 8 
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        M=varargin{4};
        Sx=varargin{5};
        Sy=varargin{6};
        Sz=varargin{7};
        c=varargin{8};
        interpMethod='linear';
    case 9
        X=varargin{1};
        Y=varargin{2};
        Z=varargin{3};
        M=varargin{4};
        Sx=varargin{5};
        Sy=varargin{6};
        Sz=varargin{7};
        c=varargin{8};
        interpMethod=varargin{9};
end
%% Create contours
% Open figure window with visibility off

ht=figure('Visible','off');
if ndims(M)==3
    
else %2D array
    X(:,:,2)=X; Y(:,:,2)=Y; Z(:,:,2)=Z+1; M(:,:,2)=M;
end

%Contourslice is used since it provides children in handles
h=contourslice(X,Y,Z,M,Sx,Sy,Sz,[c c],interpMethod); %This will plot in invisible window

%% Parse children
% 

C=cell(1,numel(h));
for i=1:length(h)
    %Get coordinates
    x=get(h(i),'XData'); y=get(h(i),'YData'); z=get(h(i),'ZData');
    
    %Vertices matrix
    V=[x(:) y(:) z(:)];
    V=V(all(~isnan(V),2),:); %Removing NaN's
    
    %Store in cell array
    C{i}=V;
end

%Get group sizes
group_sizes = cell2mat(cellfun(@(x) size(x,1), C,'UniformOutput',0)');

%Remove contours of lenght 1
C=C(group_sizes>1);
group_sizes=group_sizes(group_sizes>1);

%Sort C array according to group size
[group_sizes,ind_sort]=sort(group_sizes);
C=C(ind_sort);

close(ht);

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
