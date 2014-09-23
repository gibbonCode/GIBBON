function [V,F]=voronoiDiagramEven2D(varargin)

%% DERIVE STANDARD 2D VORONOI DIAGRAM

%Get Delaunay tesselation
DT=varargin{1};

%Computer standard Voronoi output containing vertices matrix and faces cell
[Vv,Fv_cell] =voronoiDiagram(DT);
numCell=numel(Fv_cell);

%Replace inf entries by NaN
Vv(isinf(Vv))=NaN;

%% Creating even Voronoi cell sampling

%Get number of points 
switch nargin    
    case 2
        np=varargin{2};
        interpMethod='linear';
    case 3
        np=varargin{2};
        interpMethod=varargin{3};
    otherwise
        np=[];
        interpMethod='linear';
end

%Use maximum occuring number of points if np=[]
if isempty(np)
    np=max(cellfun(@numel,Fv_cell));
end

%Resampling Voronoi cells
V=zeros(numCell*np,2); %Vertex matrix
F=reshape([1:size(V,1)]',np,size(V,1)/np)'; %Face matrix
for q=1:1:numCell
    
    %Get current Voronoi cell vertices   
    vv=Vv(Fv_cell{q},:);
    
    %Resample cell all vertices are valid
    if any(isnan(vv))
        vv=NaN;
    else        
        [vv]=evenlySampleCurve(vv,np,interpMethod,1);
    end
    
    %Store new points in array
    V(F(q,:),:)=vv;        
    
end

end