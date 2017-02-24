function axisGeom(varargin)

switch nargin
    case 0
        h=gca;
        fontSize=15;
    case 1
        h=varargin{1};
        fontSize=15;
    case 2
        h=varargin{1};
        fontSize=varargin{2};
end

if isempty(h)
    h=gca;
end

axes(h);
view(3);
xlabel('X','FontSize',fontSize);ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
set(gca,'FontSize',fontSize);
axis equal; axis vis3d; axis tight;
grid on; box on; 
