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


