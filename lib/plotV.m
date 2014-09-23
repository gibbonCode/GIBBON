function varargout=plotV(V,varargin)

nDims=size(V,2);

%Ad zeros if input is 2D
if nDims==2;
    V(:,3)=0; 
end

hp=plot3(V(:,1),V(:,2),V(:,3),varargin{:});

if nDims==2
    view(2);
else
    view(3); 
end

% varargout{1}=hp;
switch nargout
    case 1
        varargout{1}=hp;
end
   



