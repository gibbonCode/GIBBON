function [C,C_siz]=contour_group(X,Y,Z,M,Sx,Sy,Sz,c)

%% Create contours
% Open figure window with visibility off

ht=figure('Visible','off');
if ndims(M)==3
    
else %2D array
    X(:,:,2)=X; Y(:,:,2)=Y; Z(:,:,2)=Z+1; M(:,:,2)=M;
end

%Contourslice is used since it provides children in handles
h=contourslice(X,Y,Z,M,Sx,Sy,Sz,[c c]); %This will plot in invisible window

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
group_sizes = cell2mat(cellfun(@(x) numel(x), C,'UniformOutput',0)');

%Sort C array according to group size
[C_siz,ind_sort]=sort(group_sizes);
C=C(ind_sort);

close(ht);

end


