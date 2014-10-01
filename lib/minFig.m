function minFig(varargin)

%% Parse input

switch nargin
    case 0
        hf=gcf;
    case 1
        hf=varargin{1};
end

if isempty(hf)
    hf=gcf; %Take the current figure window handle if the handle is empty
end

%% Turn off warning (JavaFrame will become obsolete, this will lead to error in the future)
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

%% Minimize using JaveFrame
try
    pause(0.01);
    jFrame=get(hf,'JavaFrame');
%     set(jFrame,'Minimized',1);
    jFrame.setMinimized(1);
catch exception    
    rethrow(exception)
end

